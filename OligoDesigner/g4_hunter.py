import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict, Optional, Union


class G4Hunter:
    def __init__(self, window_size: int = 20, score_threshold: float = 0.5):
        """
        Initialize G4Hunter with parameters.
        :param window_size (int): Length of sliding window for scoring.
        :param score_threshold (float): Absolute threshold to consider potential G4.
        """
        self.window_size = window_size
        self.score_threshold = score_threshold

    def ReadSeq(self, seq: str) -> List[str]:
        """
        Read a single DNA sequence string and return list containing it.
        :param seq (str): DNA sequence (A/T/G/C, case-insensitive).
        :return: List[str]: List containing the sequence string.
        """
        return [seq.upper()]

    def BaseScore(self, seq: str) -> Tuple[str, List[int]]:
        """
        Compute a score for each nucleotide: G/G-tracts = +1..+4, C/C-tracts = -1..-4, else 0.
        :param seq (str): DNA sequence.
        :return: Tuple[str, List[int]]: Original seq.upper, list of per-base scores.
        """
        seq = seq.upper()
        scores = []
        i = 0
        n = len(seq)

        while i < n:
            if seq[i] == 'G':
                j = i
                while j < n and seq[j] == 'G':
                    j += 1
                g_count = min(j - i, 4)
                # Assign 1..g_count to first g_count positions, 4 to remaining Gs
                for k in range(g_count):
                    scores.append(k + 1)
                for _ in range(j - i - g_count):
                    scores.append(4)
                i = j
            elif seq[i] == 'C':
                j = i
                while j < n and seq[j] == 'C':
                    j += 1
                c_count = min(j - i, 4)
                for k in range(c_count):
                    scores.append(-(k + 1))
                for _ in range(j - i - c_count):
                    scores.append(-4)
                i = j
            else:
                scores.append(0)
                i += 1

        return seq, scores

    def CalScore(self, scores: List[int], window: int) -> List[float]:
        """
        Compute sliding window average scores.
        :param scores (List[int]): Base scores from BaseScore().
        :param window (int): Window size.
        :return: List[float]: Moving average scores.
        """
        if window <= 0:
            raise ValueError("Window size must be positive.")
        if len(scores) < window:
            return []

        avg_scores = []
        for i in range(len(scores) - window + 1):
            window_scores = scores[i:i + window]
            avg_scores.append(sum(window_scores) / window)
        return avg_scores

    def GFinder(self, seq: str) -> Dict[str, Union[str, List, float]]:
        """
        Full analysis pipeline for one sequence.
        :param seq (str): Input DNA sequence.
        :return: Dict with keys:
                - 'sequence': sequence string
                - 'base_scores': list of per-base scores
                - 'window_scores': list of window-averaged scores
                - 'len': length of sequence
        """
        seq_clean, base_scores = self.BaseScore(seq)
        window_scores = self.CalScore(base_scores, self.window_size)

        return {
            "sequence": seq_clean,
            "base_scores": base_scores,
            "window_scores": window_scores,
            "length": len(seq_clean),
        }

    def GetG4(self, seq: str) -> List[Dict]:
        """
        Find putative G4 regions exceeding score threshold.
        :param seq (str): Input DNA sequence.
        :return: List[Dict]: Each dict with keys: 'start', 'end', 'length', 'score', 'sequence'
        """
        result = self.GFinder(seq)
        base_scores = result["base_scores"]
        window_scores = result["window_scores"]
        window = self.window_size
        threshold = self.score_threshold

        g4_candidates = []

        for i, score in enumerate(window_scores):
            if abs(score) >= threshold:
                start = i
                end = i + window - 1
                subseq = seq[start:end + 1].upper()
                g4_candidates.append({
                    "start": start,
                    "end": end,
                    "sequence": subseq,
                    "length": len(subseq),
                    "score": round(score, 3),
                })

        return g4_candidates

    def WriteSeq(self, seq: str) -> Dict[str, Union[str, List[Dict], float]]:
        """
        Group overlapping G4 candidates and aggregate them.
        :param seq (str): Input DNA sequence.
        :return: Dict: Key 'g4_clusters' with list of merged G4 clusters, 'sum_score' 'mean_score'.
        """
        g4s = self.GetG4(seq)
        if not g4s:
            return {"g4_clusters": [], "mean_score": 0.0, "sum_score": 0.0}

        # Merge overlapping/adjacent G4s
        clusters = []
        current = g4s[0]
        for g4 in g4s[1:]:
            if g4["start"] <= current["end"] + 1:
                current["end"] = max(current["end"], g4["end"])
                current["sequence"] = seq[current["start"]:current["end"] + 1].upper()
                current["length"] = len(current["sequence"])
                start_idx = current["start"]
                end_idx = current["end"] - self.window_size + 2
                window_scores = self.CalScore(self.BaseScore(seq)[1], self.window_size)[start_idx:end_idx + 1]
                current["score"] = round(np.mean(window_scores), 3) if window_scores else 0.0
            else:
                clusters.append(current)
                current = g4
        clusters.append(current)

        mean_score = np.mean([abs(c["score"]) for c in clusters]) if clusters else 0.0
        sum_score = np.sum([abs(c["score"]) for c in clusters]) if clusters else 0.0

        return {"g4_clusters": clusters, "mean_score": round(mean_score, 3), "sum_score": round(sum_score, 3)}

    def plot2(self, scores: List[float], show_plot: bool = True, title: str = "G4 Hunter Score Profile") -> Optional[plt.Figure]:
        """
        Plot the sliding window scores.
        :param scores (List[float]): Window scores (e.g., from GFinder()['window_scores']).
        :param show_plot (bool): Whether to display plot immediately.
        :param title (str): Plot title.
        :return: Optional[matplotlib.figure.Figure]: Figure object (useful for saving/exporting).
        """
        if not scores:
            raise ValueError("No scores to plot.")

        fig, ax = plt.subplots(figsize=(12, 4))
        ax.plot(scores, 'b-', linewidth=1)
        ax.axhline(y=self.score_threshold, color='r', linestyle='--', label=f"Threshold (+{self.score_threshold})")
        ax.axhline(y=-self.score_threshold, color='r', linestyle='--', label=f"Threshold (-{self.score_threshold})")
        ax.set_xlabel("Position (nt)")
        ax.set_ylabel("G4 Score")
        ax.set_title(title)
        ax.grid(True)
        ax.legend()

        if show_plot:
            plt.show()
        return fig

    def run(self, seq: str) -> Dict:
        """
        Convenience method: full report as dictionary.
        :param seq (str): DNA sequence.
        :return: Dict: Combined results.
        """
        gfinder_res = self.GFinder(seq)
        g4_candidates = self.GetG4(seq)
        clusters = self.WriteSeq(seq)

        return {
            "raw_sequence": gfinder_res["sequence"],
            "length": gfinder_res["length"],
            "base_scores": gfinder_res["base_scores"],
            "window_scores": gfinder_res["window_scores"],
            "putative_g4s": g4_candidates,
            "g4_clusters": clusters["g4_clusters"],
            "mean_score": clusters["mean_score"],
            "sum_score": clusters["sum_score"]
        }

if __name__ == "__main__":
    import re
    g4_pattern = re.compile(r"G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}", re.I)
    def find_g4(seq):
        return [m.span() for m in g4_pattern.finditer(seq)]


    seq1 = "GGGCGGGTGGGCGGG"
    seq2 = "acaGGGCcatctGGGTtctaGGGCcatatGGGacacaGGGCcatctGGGTtctaGGGCcatatGGGacaa"
    seq3 = "GGGTTGGGTTGGGTTGGG"
    hunter = G4Hunter(window_size=25, score_threshold=0.5)

    print(hunter.run(seq1))
    print(hunter.run(seq2)['sum_score'])
    print(hunter.run(seq3))


    print(find_g4(seq1), find_g4(seq2), find_g4(seq3))
