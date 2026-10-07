#!/usr/bin/env python3
"""Build labelled SSU source sequences that are independent of the Metaxa2 database.

Sources
  bacteria, archaea, eukaryota  SILVA SSU Ref NR99 (nuclear 18S for eukaryota); SILVA
                                entries labelled Chloroplast or Mitochondria are not used
  mitochondria, chloroplast     SSU rRNA genes annotated in NCBI RefSeq mitochondrion
                                and plastid genomes (plastid: chloroplasts only)
  negatives                     SILVA LSU Ref NR99 sequences, plus composition-preserving
                                shuffles of sampled SSU sources

Independence
  Any sequence whose accession also occurs in the Metaxa2 SSU database is dropped,
  identical sequences are collapsed, and for every candidate the identity to its
  closest Metaxa2 sequence is computed with vsearch. Sources are then drawn evenly
  from three novelty bands per origin (99% or more, 97 to 99%, below 97%), filling a
  band's shortfall from the others, so that accuracy can be reported by novelty.

Outputs in --out:
  sources.fasta     source sequences (IDs src000001, ...)
  sources.tsv       source_id, origin, subgroup, database, accession, length,
                    nearest_metaxa2, identity, band
  negatives.fasta   negative sequences (IDs neg000001, ...)
  negatives.tsv     source_id, kind, database, accession, length
"""
import argparse
import gzip
import os
import random
import re
import subprocess
import sys
from collections import defaultdict

ORIGINS = ["bacteria", "archaea", "eukaryota", "mitochondria", "chloroplast"]
MIN_LEN = {"bacteria": 1200, "archaea": 1200, "eukaryota": 1500, "mitochondria": 600, "chloroplast": 1300}
BANDS = ["99-100", "97-99", "<97"]
COMP = str.maketrans("ACGTN", "TGCAN")

MITO_SSU = re.compile(r"\b(12S|s-rRNA|small[ _-]subunit|rrnS|rns|18S|rrn18|SSU)\b", re.I)
MITO_NOT = re.compile(r"\b(16S|l-rRNA|large[ _-]subunit|rrnL|rnl|23S|26S|28S|rrn26|5S|rrn5)\b", re.I)
PLASTID_SSU = re.compile(r"\b(16S|small[ _-]subunit|rrn16|rrnS|SSU)\b", re.I)
PLASTID_NOT = re.compile(r"\b(23S|large[ _-]subunit|rrn23|4\.5S|rrn4\.5|5S|rrn5)\b", re.I)


def log(msg):
    print(f">> {msg}", file=sys.stderr, flush=True)


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta(path):
    name, chunks = None, []
    with opener(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].rstrip("\n"), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def clean(seq):
    return seq.upper().replace("U", "T")


def ambiguous_fraction(seq):
    return sum(c not in "ACGT" for c in seq) / max(1, len(seq))


def accession_root(token):
    token = token.split("|")[-1]
    return token.split(".")[0]


class Reservoir:
    """Seeded reservoir sample of up to k items."""

    def __init__(self, k, rng):
        self.k, self.rng, self.items, self.seen = k, rng, [], 0

    def add(self, item):
        self.seen += 1
        if len(self.items) < self.k:
            self.items.append(item)
        else:
            j = self.rng.randrange(self.seen)
            if j < self.k:
                self.items[j] = item


def silva_ssu(path, excluded, pools):
    skipped = defaultdict(int)
    for header, seq in read_fasta(path):
        acc_token, _, tax = header.partition(" ")
        if "Chloroplast" in tax or "Mitochondria" in tax:
            skipped["organelle label"] += 1
            continue
        origin = {"Bacteria": "bacteria", "Archaea": "archaea", "Eukaryota": "eukaryota"}.get(tax.split(";")[0])
        if origin is None:
            continue
        if accession_root(acc_token) in excluded:
            skipped["in Metaxa2 database"] += 1
            continue
        seq = clean(seq)
        if len(seq) < MIN_LEN[origin] or ambiguous_fraction(seq) > 0.005:
            skipped["too short or ambiguous"] += 1
            continue
        pools[origin].add({"origin": origin, "subgroup": tax.split(";")[1] if ";" in tax else "",
                           "database": "SILVA SSU", "accession": acc_token, "seq": seq})
    log("SILVA SSU skipped: " + ", ".join(f"{k} {v}" for k, v in skipped.items()))


def parse_location(loc):
    """Return [(start, end, strand)] (0-based, end exclusive), or None for partial or
    split features (split rRNA genes are rare and their part order is ambiguous)."""
    if "<" in loc or ">" in loc or "join" in loc or "order" in loc:
        return None
    strand = -1 if loc.startswith("complement(") else 1
    parts = []
    for a, b in re.findall(r"(\d+)\.\.(\d+)", loc):
        parts.append((int(a) - 1, int(b)))
    if not parts:
        return None
    return [(a, b, strand) for a, b in parts]


def extract(seq, parts):
    pieces = [seq[a:b] for a, b, _ in parts]
    s = "".join(pieces)
    if parts[0][2] == -1:
        s = s.translate(COMP)[::-1]
    return s


def genbank_records(path):
    """Yield (accession, lineage, organelle, [(location, qualifiers)], sequence) per record."""
    rec = None
    with opener(path) as f:
        for line in f:
            if line.startswith("LOCUS"):
                rec = {"acc": "", "lineage": [], "organelle": "", "features": [], "seq": [],
                       "section": None}
            elif rec is None:
                continue
            elif line.startswith("ACCESSION"):
                rec["acc"] = line.split()[1]
            elif line.startswith("  ORGANISM"):
                rec["section"] = "organism"
            elif rec["section"] == "organism" and line.startswith("            "):
                rec["lineage"].append(line.strip())
            elif line.startswith("FEATURES"):
                rec["section"] = "features"
            elif line.startswith("ORIGIN"):
                rec["section"] = "origin"
            elif line.startswith("//"):
                yield (rec["acc"], " ".join(rec["lineage"]), rec["organelle"], rec["features"],
                       "".join(rec["seq"]).upper())
                rec = None
            elif rec["section"] == "features" and line[5:6] != " " and len(line) > 21:
                key, loc = line[5:21].strip(), line[21:].strip()
                rec["features"].append([key, loc, {}, None])
            elif rec["section"] == "features" and line.startswith("                     "):
                text = line[21:].rstrip("\n")
                feat = rec["features"][-1] if rec["features"] else None
                if feat is None:
                    continue
                if text.startswith("/"):
                    k, _, v = text[1:].partition("=")
                    feat[2][k] = v.strip('"')
                    feat[3] = k
                elif feat[3] is None:
                    feat[1] += text.strip()        # location continues
                else:
                    feat[2][feat[3]] += " " + text.strip().strip('"')
                if feat[0] == "source" and "organelle" in feat[2]:
                    rec["organelle"] = feat[2]["organelle"]
            elif rec["section"] == "origin":
                rec["seq"].append(re.sub(r"[^A-Za-z]", "", line))
            elif rec["section"] == "organism" and not line.startswith("            "):
                rec["section"] = None


def refseq_organelles(paths, excluded, pools):
    counts = defaultdict(int)
    for path in paths:
        for acc, lineage, organelle, features, seq in genbank_records(path):
            org = organelle.lower()
            if "mitochondrion" in org:
                origin, ok, bad = "mitochondria", MITO_SSU, MITO_NOT
                subgroup = "metazoan" if "Metazoa" in lineage else "non-metazoan"
            elif "chloroplast" in org:
                origin, ok, bad = "chloroplast", PLASTID_SSU, PLASTID_NOT
                subgroup = "chloroplast"
            else:
                continue
            if accession_root(acc) in excluded:
                counts["genome in Metaxa2 database"] += 1
                continue
            for key, loc, quals, _ in features:
                if key != "rRNA":
                    continue
                label = " ".join([quals.get("product", ""), quals.get("gene", "")])
                if not ok.search(label) or bad.search(label):
                    continue
                parts = parse_location(loc)
                if parts is None:
                    counts["partial or split feature"] += 1
                    continue
                s = clean(extract(seq, parts))
                if len(s) < MIN_LEN[origin] or ambiguous_fraction(s) > 0.005:
                    counts["too short or ambiguous"] += 1
                    continue
                counts[f"{origin} SSU genes"] += 1
                pools[origin].add({"origin": origin, "subgroup": subgroup, "database": "RefSeq",
                                   "accession": f"{acc}:{loc}", "seq": s})
    log("RefSeq organelles: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))


def nearest_identity(cands, metaxa2_fasta, workdir, threads, label):
    query = os.path.join(workdir, f"pool_{label}.fasta")
    with open(query, "w") as f:
        for i, c in enumerate(cands):
            f.write(f">c{i}\n{c['seq']}\n")
    hits = os.path.join(workdir, f"pool_{label}_vs_metaxa2.tsv")
    subprocess.run(["vsearch", "--usearch_global", query, "--db", metaxa2_fasta, "--id", "0.5",
                    "--strand", "both", "--maxaccepts", "10", "--maxrejects", "64",
                    "--threads", str(threads), "--blast6out", hits, "--output_no_hits",
                    "--quiet"], check=True)
    best = {}
    with open(hits) as f:
        for line in f:
            q, t, pid = line.split("\t")[:3]
            pid = float(pid)
            if t != "*" and pid > best.get(q, ("", -1.0))[1]:
                best[q] = (t, pid)
    for i, c in enumerate(cands):
        t, pid = best.get(f"c{i}", ("", 0.0))
        c["nearest"], c["identity"] = t, pid
        c["band"] = "99-100" if pid >= 99 else "97-99" if pid >= 97 else "<97"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--silva-ssu", required=True)
    ap.add_argument("--silva-lsu", required=True)
    ap.add_argument("--refseq", nargs="+", required=True, help="RefSeq mitochondrion and plastid .gbff(.gz) files")
    ap.add_argument("--metaxa2-fasta", required=True, help="Metaxa2 SSU database sequences as FASTA")
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-origin", type=int, default=2000)
    ap.add_argument("--pool-factor", type=int, default=5)
    ap.add_argument("--negatives", type=int, default=2000)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)

    excluded = {accession_root(h.split()[0]) for h, _ in read_fasta(a.metaxa2_fasta)}
    log(f"{len(excluded)} accessions in the Metaxa2 SSU database will be excluded")

    pool_k = a.per_origin * a.pool_factor
    pools = {o: Reservoir(pool_k, random.Random(f"{a.seed}-{o}")) for o in ORIGINS}
    silva_ssu(a.silva_ssu, excluded, pools)
    refseq_organelles(a.refseq, excluded, pools)

    chosen = []
    seen_seqs = set()
    for o in ORIGINS:
        cands = []
        for c in pools[o].items:
            if c["seq"] not in seen_seqs:
                seen_seqs.add(c["seq"])
                cands.append(c)
        log(f"{o}: {pools[o].seen} candidates, pool of {len(cands)} unique sequences")
        if not cands:
            continue
        nearest_identity(cands, a.metaxa2_fasta, a.out, a.threads, o)
        by_band = {b: [c for c in cands if c["band"] == b] for b in BANDS}
        for b in BANDS:
            rng.shuffle(by_band[b])
        target = a.per_origin // len(BANDS)
        take = {b: min(target, len(by_band[b])) for b in BANDS}
        spare = a.per_origin - sum(take.values())
        for b in BANDS:                       # fill shortfalls from bands that have more
            extra = min(spare, len(by_band[b]) - take[b])
            take[b] += extra
            spare -= extra
        picked = [c for b in BANDS for c in by_band[b][:take[b]]]
        log(f"{o}: picked {len(picked)} ("
            + ", ".join(f"{b} {take[b]} of {len(by_band[b])}" for b in BANDS) + ")")
        chosen.extend(picked)

    with open(os.path.join(a.out, "sources.fasta"), "w") as fa, \
         open(os.path.join(a.out, "sources.tsv"), "w") as tsv:
        tsv.write("source_id\torigin\tsubgroup\tdatabase\taccession\tlength\tnearest_metaxa2\tidentity\tband\n")
        for i, c in enumerate(chosen, 1):
            sid = f"src{i:06d}"
            fa.write(f">{sid}\n{c['seq']}\n")
            tsv.write(f"{sid}\t{c['origin']}\t{c['subgroup']}\t{c['database']}\t{c['accession']}\t"
                      f"{len(c['seq'])}\t{c['nearest']}\t{c['identity']:.1f}\t{c['band']}\n")

    # negatives: half LSU, half shuffled SSU sources
    n_lsu = a.negatives // 2
    lsu = Reservoir(n_lsu, random.Random(f"{a.seed}-lsu"))
    for header, seq in read_fasta(a.silva_lsu):
        seq = clean(seq)
        if len(seq) >= 1500 and ambiguous_fraction(seq) <= 0.005:
            lsu.add((header.split()[0], seq))
    shuffled = []
    for c in rng.sample(chosen, min(a.negatives - len(lsu.items), len(chosen))):
        s = list(c["seq"])
        rng.shuffle(s)
        shuffled.append((c["accession"], "".join(s)))
    with open(os.path.join(a.out, "negatives.fasta"), "w") as fa, \
         open(os.path.join(a.out, "negatives.tsv"), "w") as tsv:
        tsv.write("source_id\tkind\tdatabase\taccession\tlength\n")
        n = 0
        for kind, db, items in (("LSU", "SILVA LSU", lsu.items), ("shuffled SSU", "shuffled", shuffled)):
            for acc, seq in items:
                n += 1
                fa.write(f">neg{n:06d}\n{seq}\n")
                tsv.write(f"neg{n:06d}\t{kind}\t{db}\t{acc}\t{len(seq)}\n")
    log(f"wrote {len(chosen)} sources and {n} negatives to {a.out}")


if __name__ == "__main__":
    main()
