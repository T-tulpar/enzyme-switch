#!/usr/bin/env python3
"""codon_optimization.py — codon optimization + CAI analysis (Module 4).


INPUT  : data/petase_seq.fasta          (VERIFIED PETase sequence, you fetch)
         data/ecoli_k12_codon_usage.txt (E. coli K-12 codon usage, you fetch)
OUTPUT : data/petase_cds_optimized.fa + figures/cai_comparison.png

CAI = geometric mean of Sharp & Li-style weights = count / max count
within each synonymous group.
"""

import math
import pathlib
import re

import matplotlib
matplotlib.use("Agg")        # non-interactive backend: no display required
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "data"
FIG_DIR = REPO_ROOT / "figures"
FIG_DIR.mkdir(exist_ok=True)

INPUT_FASTA = DATA_DIR / "petase_seq.fasta"
CODON_TABLE = DATA_DIR / "ecoli_k12_codon_usage.txt"
OUTPUT_FASTA = DATA_DIR / "petase_cds_optimized.fa"

RESTRICTION_SITES = ("CCATGG", "CTCGAG")  # NcoI, XhoI

# synonymous codon groups, bacterial codon table (NCBI #11)
TABLE = CodonTable.unambiguous_dna_by_id[11]
SYN = {}
for codon, aa in TABLE.forward_table.items():
    SYN.setdefault(aa, []).append(codon)


def gc_content(seq):
    """GC percentage, computed by hand (version-proof)."""
    seq = str(seq).upper()
    return (seq.count("G") + seq.count("C")) / len(seq) * 100.0


def load_codon_usage(path):
    """Parse Kazusa-style table into {codon: relative frequency}.

    Tolerant regex: grabs 3-letter codons (A/C/G/U) + the number that
    follows (frequency per thousand OR count — either works because we
    only use relative values within each group). U -> T (RNA -> DNA).
    """
    text = path.read_text()
    items = re.findall(r"([ACGU]{3})\s+(?:[A-Z]\s+)?([\d.]+)", text)
    usage = {c.replace("U", "T"): float(n) for c, n in items}
    if len(usage) < 60:
        raise ValueError(
            f"Parsed only {len(usage)} codons from {path.name} "
            "(want 64) - file may not be in Kazusa format."
        )
    return usage


def build_weights(usage):
    """Sharp & Li-style weight = count / max count within the group."""
    w = {}
    for aa, codons in SYN.items():
        vals = [usage.get(c, 0.0) for c in codons]
        m = max(vals)
        for c, v in zip(codons, vals):
            w[c] = (v / m) if m > 0 else 0.0
    return w


def cai_of(cds, weights):
    """Geometric mean of codon weights over the CDS."""
    codons = [cds[i:i+3] for i in range(0, len(cds) - 3, 3)]
    ws = [weights[c] for c in codons if weights.get(c, 0.0) > 0]
    if not ws:
        return float("nan")
    return math.exp(sum(math.log(x) for x in ws) / len(ws))


def introduces_site(prefix, codon, sites):
    """True if appending this codon would create a restriction site that
    spans the boundary. (A site fully inside the existing prefix already
    existed, so we only need to inspect the last 8 characters.)"""
    window = (prefix + codon)[-8:].upper()
    return any(s in window for s in sites)

def optimize_protein(protein, usage):
    """Pick the best synonymous codon per amino acid, avoiding internal
    restriction sites, using BACKTRACKING so a dead-end rolls back and
    tries a different earlier codon instead of giving up."""
    sites = RESTRICTION_SITES
    # candidate codons per position, best (highest usage) first
    candidates = [
        sorted(SYN.get(aa, []), key=lambda c: -usage.get(c, 0.0))
        for aa in protein
    ]
    chosen = []                       # codons chosen so far

    def solve(i):
        if i == len(protein):         # reached the end -> success
            return True
        for c in candidates[i]:       # try each codon, best first
            if introduces_site("".join(chosen), c, sites):
                continue              # would create a cut site -> skip
            chosen.append(c)
            if solve(i + 1):          # does the rest work out?
                return True
            chosen.pop()              # didn't -> undo and try the next
        return False                  # nothing worked here -> backtrack further

    if not solve(0):
        raise ValueError(
            "Could not assemble a CDS avoiding NcoI/XhoI. "
            "Relax the constraint or hand-design the problem region."
        )
    return "".join(chosen)

def load_input():
    if not INPUT_FASTA.exists():
        print("[STOP] Input file not found:", INPUT_FASTA)
        print("  Fetch the VERIFIED PETase sequence first (integrity policy).")
        print("  Expected source: UniProt A0A0K8P6T7 - VERIFY at fetch!")
        return None, None
    record = SeqIO.read(INPUT_FASTA, "fasta")
    seq = str(record.seq).upper()
    if set(seq) <= set("ACGT"):              # looks like nucleotide CDS
        dna = seq
        protein = str(Seq(dna).translate())
        if protein.endswith("*"):
            protein = protein[:-1]
        return dna, protein
    return None, seq                          # assume protein input


def main():
    dna, protein = load_input()
    if protein is None:
        return

    if not CODON_TABLE.exists():
        print("[STOP] Codon table not found:", CODON_TABLE)
        print("  Download the E. coli K-12 codon-usage table from Kazusa")
        print("  and save it as this file.")
        return
    usage = load_codon_usage(CODON_TABLE)
    weights = build_weights(usage)

    optimized = optimize_protein(protein, usage)

    # ---- THE critical check: protein must be identical -----------------
    back = str(Seq(optimized).translate())
    if back.endswith("*"):
        back = back[:-1]
    ok = back == protein
    print("=== Back-translation validation ===")
    print(f"  protein length : {len(protein)} aa")
    print(f"  optimized CDS  : {len(optimized)} bp (frame OK: {len(optimized) % 3 == 0})")
    print(f"  identical?     : {'PASS - no mutations introduced' if ok else 'FAIL!'}")
    if not ok:
        raise SystemExit(1)

    cai_after = cai_of(optimized, weights)
    print("\n=== CAI / GC summary (real E. coli K-12 table) ===")
    print(f"  Optimized CDS : CAI = {cai_after:.3f} | GC = {gc_content(optimized):5.1f}%")
    if dna is not None:
        cai_before = cai_of(dna, weights)
        print(f"  Original CDS  : CAI = {cai_before:.3f} | GC = {gc_content(dna):5.1f}%")
        print(f"  CAI gain      : {100 * (cai_after - cai_before) / cai_before:+.1f}%")

    hits = [s for s in RESTRICTION_SITES if s in optimized]
    print(f"  Internal {RESTRICTION_SITES} sites: {'NONE (good)' if not hits else hits}")

    rec = SeqRecord(Seq(optimized), id="PETase_Ec_opt",
                    description="codon-optimized for E. coli K-12 (in silico)")
    SeqIO.write(rec, OUTPUT_FASTA, "fasta")
    print(f"  Saved: {OUTPUT_FASTA}")

    # figure
    fig, ax = plt.subplots(figsize=(6, 4))
    labels, values = ["Optimized"], [cai_after]
    if dna is not None:
        labels, values = ["Original", "Optimized"], [cai_of(dna, weights), cai_after]
    ax.bar(labels, values, color=["gray", "crimson"])
    ax.set_ylabel("CAI (real E. coli K-12 table)")
    ax.set_ylim(0, 1.05)
    ax.set_title("Codon Adaptation Index (in silico)")
    fig.tight_layout()
    out = FIG_DIR / "cai_comparison.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved figure: {out}")


if __name__ == "__main__":
    main()
