"""Count influenza (LAIV) reads in nasal RNA-seq, to derive a replication endpoint.

Pre-check and endpoint for H1/H4 of docs/PREREG_vaccine_axes.md when the published
PCR shedding data cannot be linked to the transcriptomes. Vaccine replication is
instead measured in the authors' own deposited reads, matched to participants by
their own sample titles within GSE230494.

Exact k-mer screen (no aligner): a read counts as influenza if any of a few
non-overlapping k-mer windows matches a k-mer of the reference set (either
strand). With k=25 a chance match against human sequence is ~1e-15 per lookup, so
positives are specific; sensitivity is limited by strain divergence and
sequencing error, which only ever makes the count conservative.

Usage: python validation/laiv_viral_reads.py <refs.fa> <fastq_dir> <out.csv> [--max-reads N]
"""
from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

K = 25
STEP = 25                      # non-overlapping windows within a 75 bp read
COMP = bytes.maketrans(b"ACGTN", b"TGCAN")


def read_fasta(path):
    name, buf = None, []
    for line in open(path, "rb"):
        if line.startswith(b">"):
            if name:
                yield name, b"".join(buf)
            name, buf = line[1:].strip().decode("utf8", "replace"), []
        else:
            buf.append(line.strip().upper())
    if name:
        yield name, b"".join(buf)


def build_kmers(path):
    ks = set()
    n_seq = 0
    for _, seq in read_fasta(path):
        n_seq += 1
        rc = seq.translate(COMP)[::-1]
        for s in (seq, rc):
            for i in range(len(s) - K + 1):
                km = s[i:i + K]
                if b"N" not in km:
                    ks.add(km)
    return ks, n_seq


def screen(path, kmers, max_reads=None):
    total = hits = 0
    with gzip.open(path, "rb") as fh:
        try:
            for i, line in enumerate(fh):
                if i & 3 != 1:                      # sequence is every 4th line, offset 1
                    continue
                total += 1
                s = line.strip()
                for j in range(0, len(s) - K + 1, STEP):
                    if s[j:j + K] in kmers:
                        hits += 1
                        break
                if max_reads and total >= max_reads:
                    break
        except (EOFError, OSError) as exc:          # truncated .gz (partial download)
            print(f"  [truncated stream after {total} reads: {exc}]", file=sys.stderr)
    return total, hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("refs"); ap.add_argument("fastq_dir"); ap.add_argument("out")
    ap.add_argument("--max-reads", type=int, default=None)
    a = ap.parse_args()

    t0 = time.time()
    kmers, n_seq = build_kmers(a.refs)
    print(f"reference: {n_seq} sequences, {len(kmers)} k-mers (k={K}), "
          f"{time.time() - t0:.0f}s", flush=True)

    rows = []
    for fq in sorted(Path(a.fastq_dir).glob("*.fastq.gz")):
        t = time.time()
        total, hits = screen(fq, kmers, a.max_reads)
        rows.append(dict(run=fq.name.split(".")[0], reads_screened=total, influenza_reads=hits,
                         per_million=hits / total * 1e6 if total else float("nan")))
        print(f"{rows[-1]['run']}: {hits} / {total} reads "
              f"({rows[-1]['per_million']:.2f} per million)  [{time.time() - t:.0f}s]", flush=True)

    import pandas as pd
    pd.DataFrame(rows).to_csv(a.out, index=False)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
