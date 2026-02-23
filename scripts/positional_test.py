import math
from collections import Counter
from itertools import product
from scipy.stats import binom
import sys

# --- IUPAC mapping ---
IUPAC = {
    "A": ["A"], "C": ["C"], "G": ["G"], "T": ["T"],
    "R": ["A", "G"], "Y": ["C", "T"],
    "S": ["G", "C"], "W": ["A", "T"],
    "K": ["G", "T"], "M": ["A", "C"],
    "B": ["C", "G", "T"], "D": ["A", "G", "T"],
    "H": ["A", "C", "T"], "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"]
}


# ---------------- FASTA READER ----------------
def read_fasta(path):
    header = None
    seq_chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq_chunks)
                header = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.upper())
        if header:
            yield header, "".join(seq_chunks)


# ---------------- MOTIF EXPANSION ----------------
def expand_motif(motif):
    pools = [IUPAC[b] for b in motif]
    return ["".join(p) for p in product(*pools)]


# ---------------- M2 FITTING ----------------
def fit_markov_m2(seq):
    """
    Fit second-order Markov model:
    - Initial distribution: dinucleotide frequencies
    - Transition: P(z | xy)
    """

    seq = seq.upper()
    L = len(seq)

    if L < 3:
        raise ValueError("Sequence too short for M2 model.")

    di = Counter(seq[i:i+2] for i in range(L - 1))
    tri = Counter(seq[i:i+3] for i in range(L - 2))

    # P(b1 b2)
    total_di = sum(di.values())
    P_init = {xy: (di[xy] + 1) / (total_di + 16)  # Laplace smoothing
              for xy in (a + b for a in "ACGT" for b in "ACGT")}

    # P(z | xy)
    P_trans = {}
    for x in "ACGT":
        for y in "ACGT":
            xy = x + y
            total_xy = sum(tri[xy + z] for z in "ACGT") + 4
            for z in "ACGT":
                P_trans[(xy, z)] = (tri[xy + z] + 1) / total_xy

    return P_init, P_trans


def prob_word_m2(word, P_init, P_trans):
    """
    Probability of word under M2.
    """

    p = P_init[word[0:2]]

    for i in range(2, len(word)):
        xy = word[i-2:i]
        z = word[i]
        p *= P_trans[(xy, z)]

    return p


def expected_count_m2(seq_len, k, words, P_init, P_trans):
    N = seq_len - k + 1
    if N <= 0:
        return 0.0

    total_p = sum(prob_word_m2(w, P_init, P_trans) for w in words)
    return N * total_p


# ---------------- COUNTING ----------------
def count_occurrences(seq, words):
    k = len(words[0])
    word_set = set(words)
    return sum(
        1 for i in range(len(seq) - k + 1)
        if seq[i:i+k] in word_set
    )


# ---------------- ANALYSIS PER SEQUENCE ----------------
def analyse_sequence_m2(seq, words, X):
    k = len(words[0])

    A = seq[:X]
    B = seq[X:]

    P_init_A, P_trans_A = fit_markov_m2(A)
    P_init_B, P_trans_B = fit_markov_m2(B)

    E_A = expected_count_m2(len(A), k, words, P_init_A, P_trans_A)
    E_B = expected_count_m2(len(B), k, words, P_init_B, P_trans_B)

    O_A = count_occurrences(A, words)
    O_B = count_occurrences(B, words)
    O_total = O_A + O_B

    if (E_A + E_B) > 0 and O_total > 0:
        p = E_A / (E_A + E_B)
        pval = binom.cdf(O_A, O_total, p)
    else:
        p = float("nan")
        pval = float("nan")

    return E_A, E_B, O_A, O_B, p, pval


# ---------------- FASTA LEVEL ----------------
def analyse_fasta_m2(fasta_path, motifs, X=10000):
    results = []
    for motif in motifs:
        words = expand_motif(motif)
    
        for header, seq in read_fasta(fasta_path):
            if len(seq) < 3:
                continue
            alpha=0.5 # to transform for log
            E_A, E_B, O_A, O_B, p, pval = analyse_sequence_m2(seq, words, X)
            delta_L = math.log((O_A + alpha)/(E_A + alpha)) - math.log((O_B + alpha)/(E_B + alpha))
            results.append({
            "id": header,
            "length": len(seq),
            "E_A": E_A,
            "E_B": E_B,
            "O_A": O_A,
            "O_B": O_B,
            "p_expected_A": p,
            "p_value_underrep_A": pval,
            "delta_L":delta_L
        })

    return results

results = analyse_fasta_m2(sys.argv[1],["GTATCC"], X=15000)

for r in results[:5]:
    print(r)

