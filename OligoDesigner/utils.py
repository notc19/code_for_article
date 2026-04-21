"""
GTF annotation parser for TS capture probe design.
Supports:
  - GENCODE v20+  (GRCh38, auto-detected or caller-supplied)
  - Ensembl v76+  (GRCh38, auto-detected or caller-supplied)

Version-aware behaviour
-----------------------
GENCODE
  v20-25  : gene_status / transcript_status present (ignored here)
  v20-35  : lincRNA, antisense, sense_intronic … used as separate biotypes
  v36+    : all lncRNA subtypes unified under "lncRNA"
  v20-37  : polymorphic_pseudogene transcript biotype
  v38+    : polymorphic_pseudogene -> protein_coding_LoF
  v20-43  : tag "basic"
  v44+    : tag "basic" renamed to "gencode_basic"
  v30+    : MANE_Select tag present
  v38+    : MANE_Plus_Clinical tag present

Ensembl
  e76-e101 : lincRNA / antisense … separate biotypes
  e102+    : lncRNA unified biotype
  e76-e104 : polymorphic_pseudogene
  e105+    : protein_coding_LoF
  e76-e111 : tag "basic"
  e112+    : tag "basic" renamed to "gencode_basic"
  e100+    : MANE_Select tag present
  e105+    : MANE_Plus_Clinical tag present
"""

from __future__ import annotations
import re
from logging import Logger

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')

def _parse_attributes(attribute: str) -> dict[str, list[str]]:
    """
    Parse the 9th GTF column into a dict {key: [value, ...]}  
    (multi-valued keys such as 'tag' accumulate all values in the list).
    Works for both Ensembl and GENCODE attribute formatting.
    :param attribute: string value of 9th GTF column
    :return: dict format of 9th GTF column
    """
    result: dict[str, list[str]] = {}
    for key, val in _ATTR_RE.findall(attribute):
        result.setdefault(key, []).append(val)
    return result


def _get_attr(attrs: dict[str, list[str]], *keys: str, default: str = "") -> str:
    """
    Return the first value for the first key that exists, or *default*
    :param attrs: dictionary of attributes
    """
    for k in keys:
        if k in attrs:
            return attrs[k][0]
    return default


def _has_tag(attrs: dict[str, list[str]], tag: str) -> bool:
    """
    Return True if 'tag' attribute list contains *tag*
    :param attrs: dictionary of attributes
    """
    return tag in attrs.get("tag", [])


def _detect_source(annot: list[str]) -> str:
    """
    Auto-detect annotation source from the first 200 annot lines.
    Returns "ensembl" or "gencode"
    :param annot: GTF annotation
    """
    for line in annot[:200]:
        col0 = line.split("\t")[0]
        if col0.startswith("chr"):
            return "gencode"
    return "ensembl"


def _biotype_key(source: str, level: str) -> str:
    """
    Return the correct attribute key for biotype depending on source for levels "transcript" or "gene"
    :param source: "gencode" or "ensembl"
    :param level: "transcript" or "gene"
    """
    if source == "gencode":
        return "transcript_type" if level == "transcript" else "gene_type"
    return "transcript_biotype" if level == "transcript" else "gene_biotype"


def _transcript_id_versioned(attrs: dict[str, list[str]], source: str) -> str:
    """
    Build a versioned transcript ID  ENSTnnn.N  regardless of source.
    - GENCODE: transcript_id already contains the version.
    - Ensembl: transcript_id + '.' + transcript_version.
    :param attrs: dictionary of attributes
    :param source: "gencode" or "ensembl"

    """
    tid = _get_attr(attrs, "transcript_id")
    if source == "ensembl" and "." not in tid:
        ver = _get_attr(attrs, "transcript_version")
        if ver:
            tid = f"{tid}.{ver}"
    return tid


def _is_utr_feature(feature: str) -> bool:
    """
    Return True for any UTR feature name used by source
    :param feature: "five_prime_utr", "three_prime_utr", or "UTR"
    """
    return feature in {"five_prime_utr", "three_prime_utr", "UTR"}


def _protein_coding_biotypes(source: str, version: int) -> frozenset[str]:
    """
    Return the set of transcript biotype values that should be treated as
    protein-coding for exon collection, accounting for version-specific renames.

    GENCODE v20-37 uses "polymorphic_pseudogene" for transcripts that are
    protein-coding in some haplotypes; v38+ renames this to
    "protein_coding_LoF".  Ensembl follows the same boundary at e105.
    Both sources always include "protein_coding" itself.

    :param source: "gencode" or "ensembl"
    :param version: integer release number (e.g. 44 for GENCODE v44, 112 for Ensembl e112)
    :return: frozenset of biotype strings considered protein-coding
    """
    base = {"protein_coding"}
    if source == "gencode":
        if version < 38:
            base.add("polymorphic_pseudogene")
        else:
            base.add("protein_coding_LoF")
    else:  # ensembl
        if version < 105:
            base.add("polymorphic_pseudogene")
        else:
            base.add("protein_coding_LoF")
    return frozenset(base)


def _basic_tag(source: str, version: int) -> str:
    """
    Return the correct GTF tag string for the "basic" transcript subset,
    accounting for the rename that occurred at GENCODE v44 / Ensembl e112.

    Prior to those releases the tag was "basic"; from those releases onward
    it became "gencode_basic".  Tags in archive files were NOT back-filled,
    so the correct string depends strictly on the version being parsed.

    :param source: "gencode" or "ensembl"
    :param version: integer release number
    :return: "gencode_basic" for GENCODE v44+ / Ensembl e112+, else "basic"
    """
    if source == "gencode":
        return "gencode_basic" if version >= 44 else "basic"
    return "gencode_basic" if version >= 112 else "basic"


def prepare_annot_for_muts(
    input_muts_list: list[str],
    annot: list[str],
    logger: Logger,
) -> dict | None:
    """
    Parse a GTF annotation and return a dict required for mutation coordinate lookup (HGVS c./n. notation).
    Supports Ensembl v110+ and GENCODE v49 GTF files (auto-detected).

    :param input_muts_list: List of mutation strings. Only parsed when at least one entry contains HGVS c. or n. notation (e.g. "GENE:c.123G>C")
    :param annot: GTF file lines, pre-stripped, comments excluded
    :param logger: object of class Logger
    :return: dict  {versioned_transcript_id: [chrom, exon_list, start_codon_coord, strand]} or None if no HGVS mutations are present in *input_muts_list*
    """
    # Only proceed when at least one mutation uses HGVS c./n. notation
    needs_parse = any(
        len(m.split(":")) > 1 and ("c." in m.split(":")[1] or "n." in m.split(":")[1])
        for m in input_muts_list
    )
    if not needs_parse:
        return None

    source = _detect_source(annot)
    transcript_dict: dict = {}

    i = 0
    while i < len(annot):
        fields = annot[i].split("\t")
        if len(fields) < 9:
            i += 1
            continue

        feature = fields[2]

        if feature == "transcript":
            attrs = _parse_attributes(fields[8])
            transcript_id = _transcript_id_versioned(attrs, source)
            transcript_chr = fields[0]
            strand = fields[6]

            exon_list: list[list[str]] = []
            start_codon_coord: str | int = 0

            k = i + 1
            while k < len(annot):
                sub_fields = annot[k].split("\t")
                if len(sub_fields) < 9:
                    k += 1
                    continue
                sub_feature = sub_fields[2]
                if sub_feature == "transcript":
                    break
                if sub_feature == "exon":
                    exon_list.append([sub_fields[3], sub_fields[4]])
                elif sub_feature == "start_codon":
                    start_codon_coord = sub_fields[3]
                k += 1

            if start_codon_coord:
                transcript_dict[transcript_id] = [
                    transcript_chr,
                    exon_list,
                    start_codon_coord,
                    strand,
                ]

            i = k  # jump past all child records already consumed
            continue

        i += 1

    logger.info("Annotation for mutations was prepared.")
    return transcript_dict


def prepare_annot_for_genes(
    annot: list[str],
    logger: Logger,
    source: str | None = None,
    version: int | None = None,
) -> tuple[dict, dict]:
    """
    Parse a GTF annotation and return dicts required for gene-region lookup.
    Exons are collected from the representative transcript, chosen by the following priority:
      1. MANE_Select transcript (preferred, present in GENCODE v30+ / Ensembl e100+).
      2. If no MANE_Select exists: any transcript carrying both
         Ensembl_canonical & gencode_primary
         (or version-appropriate basic tag — gencode_basic / basic).
      3. If neither criterion is met (e.g. very old annotation versions),
         the first protein-coding transcript encountered.        

    Protein-coding biotypes accepted for exon collection:
    +-----------+----------+----------------------------------------------+
    | Source    | Version  | Biotypes included                            |
    +===========+==========+==============================================+
    | GENCODE   | v20-v37  | protein_coding, polymorphic_pseudogene       |
    | GENCODE   | v38+     | protein_coding, protein_coding_LoF           |
    | Ensembl   | e76-e104 | protein_coding, polymorphic_pseudogene       |
    | Ensembl   | e105+    | protein_coding, protein_coding_LoF           |
    +-----------+----------+----------------------------------------------+

    The "basic" tag name changed at GENCODE v44 / Ensembl e112 from "basic" to "gencode_basic"; this function uses the correct name automatically
    the correct string is resolved automatically via _basic_tag(source, version).

    Returned structures
    -------------------
    gene_dict : dict
        {gene_id: [chrom, gene_coords, exon_list, utr_list,
                     exon_id_list, mane_select_transcripts,
                     mane_plus_clinical_transcripts]}

        - *gene_id*                   – unversioned Ensembl stable ID (e.g. "ENSG00000160087").
        - *chrom*                     – chromosome / scaffold name exactly as it appears in the GTF ("chr1" for GENCODE, "1" for Ensembl).
        - *gene_coords*               – [start, end] strings from the gene feature line (1-based, inclusive).
        - *exon_list*                 – [[start, end], …] ints, sorted ascending, from the selected transcript only,
                                        with stop-codon boundary trimming applied.
        - *utr_list*                  – [[start, end], …] strings, deduplicated by coordinate pair.
        - *exon_id_list*              – Ensembl exon stable IDs ("ENSE…") parallel to exon_list.
        - *mane_select_transcripts*   – versioned transcript IDs tagged MANE_Select; empty list for versions that pre-date MANE
                                        (GENCODE <v30 / Ensembl <e100).
        - *mane_plus_clinical_transcripts* – versioned transcript IDs tagged MANE_Plus_Clinical; empty list for versions that pre-date the tag
                                        (GENCODE <v38 / Ensembl <e105).
        - *transcript_list*           – all versioned transcript IDs 
        - *selected_transcript_id*         – the versioned transcript ID whose exons populate exon_list
                                             (useful for debugging / downstream reporting).                    

    gene_name_dict : dict
        {gene_symbol: gene_id} mapping HGNC-approved gene symbol to the
        unversioned Ensembl stable gene ID.  Falls back to *gene_id* as the
        symbol when gene_name is absent from the attributes.

    -----------------------
    :param annot: GTF file lines, pre-stripped, comment lines excluded
    :param logger: object of class Logger
    :param source: annotation source "gencode" or "ensembl"
                   When *None*, auto-detected from chromosome name formatting
                   ("chr" prefix → GENCODE)
    :param version: integer release number  
                    When *None*, defaults to 49 (GENCODE) or 115 (Ensembl) after source detection
    :return: tuple (gene_dict, gene_name_dict)
    """
    if source is None:
        source = _detect_source(annot)

    if version is None:
        version = 49 if source == "gencode" else 115

    tx_biotype_key = _biotype_key(source, "transcript")
    pc_biotypes = _protein_coding_biotypes(source, version)
    basic_tag = _basic_tag(source, version)

    # Determine whether MANE tags can be present in this version
    mane_select_available = (
        (source == "gencode" and version >= 30)
        or (source == "ensembl" and version >= 100)
    )
    mane_plus_available = (
        (source == "gencode" and version >= 38)
        or (source == "ensembl" and version >= 105)
    )
    # Ensembl_canonical introduced ~e109 / GENCODE v41
    canonical_available = (
        (source == "gencode" and version >= 41)
        or (source == "ensembl" and version >= 109)
    )
    # gencode_primary introduced at GENCODE v44 / Ensembl e113
    primary_available = (
        (source == "gencode" and version >= 44)
        or (source == "ensembl" and version >= 113)
    )

    gene_dict: dict = {}
    gene_name_dict: dict = {}

    i = 0
    while i < len(annot):
        fields = annot[i].split("\t")
        if len(fields) < 9:
            i += 1
            continue

        feature = fields[2]

        if feature != "gene":
            i += 1
            continue

        gene_attrs = _parse_attributes(fields[8])
        gene_id_raw = _get_attr(gene_attrs, "gene_id")
        # Strip version suffix so the key is always a stable unversioned ID
        gene_id = gene_id_raw.split(".")[0]
        gene_name = _get_attr(gene_attrs, "gene_name", default=gene_id)
        gene_chr = fields[0]
        gene_coords = [fields[3], fields[4]]

        gene_name_dict[gene_name] = gene_id

        current_tid: str | None = None
        transcript_id_list: list[list[str]] =[]
        tx_buffer: dict[str, dict] = {}
        tx_order: list[str] = []
        mane_select: list[str] = []
        mane_plus_clinical: list[str] = []

        k = i + 1
        while k < len(annot):
            sub_fields = annot[k].split("\t")
            if len(sub_fields) < 9:
                k += 1
                continue
            sub_feature = sub_fields[2]

            if sub_feature == "gene":
                break  # next gene record

            sub_attrs = _parse_attributes(sub_fields[8])

            # -- transcript-level records + MANE tag collection ---------------
            if sub_feature == "transcript":
                tid = _transcript_id_versioned(sub_attrs, source)
                current_tid = tid
                transcript_id_list.append(tid)
                
                # Collect tags for this transcript into a set for easy lookup.
                tags: set[str] = set(sub_attrs.get("tag", []))
                biotype = _get_attr(sub_attrs, tx_biotype_key)
                tx_buffer[tid] = {
                    "exons": [],        # [[start, end], ...]
                    "exon_ids": [],     # [exon_id, ...]
                    "utrs": [],         # [[start, end], ...]
                    "utr_set": set(),
                    "stop_codons": [],  # [int, ...]
                    "tags": tags,
                    "biotype": biotype,
                    "strand": sub_fields[6],
                }
                tx_order.append(tid)

                if mane_select_available and _has_tag(sub_attrs, "MANE_Select"):
                    mane_select.append(tid)
                if mane_plus_available and _has_tag(sub_attrs, "MANE_Plus_Clinical"):
                    mane_plus_clinical.append(tid)

            # -- exon: buffer into the current transcript -------------------
            elif sub_feature == "exon" and current_tid is not None:
                buf = tx_buffer[current_tid]
                biotype = _get_attr(sub_attrs, tx_biotype_key)
                if biotype in pc_biotypes:
                    exon_id = _get_attr(sub_attrs, "exon_id")
                    if exon_id and exon_id not in buf["exon_ids"]:
                        buf["exons"].append([int(sub_fields[3]), int(sub_fields[4])])
                        if source == "gencode":
                            exon_id = exon_id.split(".")[0]
                        buf["exon_ids"].append(exon_id)

            # -- UTR: buffer into the current transcript --------------------
            elif _is_utr_feature(sub_feature) and current_tid is not None:
                buf = tx_buffer[current_tid]
                utr_coords = (sub_fields[3], sub_fields[4])
                if utr_coords not in buf["utr_set"]:
                    buf["utrs"].append(list(utr_coords))
                    buf["utr_set"].add(utr_coords)

            # -- stop codon: buffer into the current transcript -------------
            elif sub_feature == "stop_codon" and current_tid is not None:
                tx_buffer[current_tid]["stop_codons"].append(int(sub_fields[4]))

            k += 1

        # ---- select the representative transcript -------------------------
        # Priority 1: MANE_Select
        selected_tid: str | None = mane_select[0] if mane_select else None

        # Priority 2: Ensembl_canonical AND (gencode_primary OR basic tag)
        if selected_tid is None:
            for tid in tx_order:
                buf = tx_buffer[tid]
                tags = buf["tags"]
                has_canonical = (not canonical_available) or ("Ensembl_canonical" in tags)
                has_primary = (
                    (primary_available and "gencode_primary" in tags)
                    or (basic_tag in tags)
                )
                if has_canonical and has_primary and buf["biotype"] in pc_biotypes:
                    selected_tid = tid
                    break

        # Priority 3: first protein-coding transcript (last-resort fallback)
        if selected_tid is None:
            for tid in tx_order:
                if tx_buffer[tid]["biotype"] in pc_biotypes:
                    selected_tid = tid
                    logger.debug(
                        f"{gene_id}: no MANE_Select, canonical or primary transcript found; "
                        f"falling back to first protein-coding transcript {tid}."
                    )
                    break

        # ---- extract exons/UTRs from the selected transcript --------------
        if selected_tid is not None and tx_buffer[selected_tid]["exons"]:
            buf = tx_buffer[selected_tid]
            sel_strand = buf["strand"]
            exon_list = buf["exons"]
            exon_id_list = buf["exon_ids"]
            utr_list = buf["utrs"]
            stop_codon_coords = buf["stop_codons"]

            # Sort exons and keep exon_id_list in the same order.
            paired = sorted(zip(exon_list, exon_id_list), key=lambda x: x[0][0])
            exon_list = [p[0] for p in paired]
            exon_id_list = [p[1] for p in paired]

            # Apply stop-codon boundary trimming.
            if stop_codon_coords:
                if sel_strand == "+":
                    for stop in stop_codon_coords:
                        for exon in exon_list:
                            if exon[0] <= stop <= exon[1]:
                                exon[1] = stop
                elif sel_strand == "-":
                    for stop in stop_codon_coords:
                        for exon in exon_list:
                            if exon[0] <= stop <= exon[1]:
                                exon[0] = stop
        else:
            exon_list = []
            exon_id_list = []
            utr_list = []
            if selected_tid is None:
                logger.debug(f"{gene_id}: no protein-coding transcript found; exon_list will be empty.")

        gene_dict[gene_id] = [
            gene_chr,
            gene_coords,
            exon_list,
            utr_list,
            exon_id_list,
            mane_select,
            mane_plus_clinical,
            transcript_id_list,
            selected_tid,
        ]

        i = k

    logger.info(
        f"Annotation (source={source}, version={version}) for genes was prepared, {len(gene_dict)} genes parsed."
    )
    return gene_dict, gene_name_dict


if __name__ == "__main__":
    import sys
    import logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    test_gtf = "../tests/FileTest_cardio_gencode_annot.gtf"

    try:
        with open(test_gtf) as fh:
            annot = [line.strip() for line in fh if not line.startswith("#")]
    except FileNotFoundError:
        print(f"Test file not found: {test_gtf}", file=sys.stderr)
        sys.exit(1)

    result = prepare_annot_for_muts(['ENST00000641515.2:c.4A>C', 'ENST00000269834.1:c.1398G>A', 'ENSG00000130164:c.123G>C'], annot, logger)

    print("annot_for_genes")
    gene_dict, gene_name_dict = prepare_annot_for_genes(annot, logger)
    #print(gene_dict['ENSG00000130164'])
    for gene_id, info in list(gene_dict.items()):
        chrom, coords, exons, utrs, exon_ids, mane_sel, mane_plus, transcripts, selected_tid = info
        print(
            f"{gene_id}  chr={chrom}  coords={coords} "
            f"exons={len(exons)} transcript_selected={selected_tid} "
            f"MANE_Select={mane_sel}\n{exon_ids}\n"
        )
