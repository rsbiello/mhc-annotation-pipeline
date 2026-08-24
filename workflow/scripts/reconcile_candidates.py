from collections import defaultdict
from pathlib import Path

blast_path = Path(snakemake.input.blast)
miniprot_path = Path(snakemake.input.miniprot)
fai_path = Path(snakemake.input.fai)

cluster_distance = int(snakemake.params.cluster_distance)
overlap_flank = int(snakemake.params.overlap_flank)

for output in snakemake.output:
    Path(output).parent.mkdir(parents=True, exist_ok=True)

chrom_lengths = {}
with fai_path.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        chrom_lengths[fields[0]] = int(fields[1])

hsps = []
with blast_path.open() as handle:
    for line in handle:
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 13:
            raise ValueError("Expected 13-column TBLASTN output; check the Snakefile outfmt")
        query, chrom = f[0], f[1]
        s1, s2 = int(f[7]), int(f[8])
        start, end = min(s1, s2) - 1, max(s1, s2)
        strand = "+" if s1 <= s2 else "-"
        family = query.split("|", 1)[0]
        hsps.append({"chrom": chrom, "start": start, "end": end, "strand": strand,
                     "query": query, "family": family, "evalue": f[10], "bitscore": f[11]})

hsps.sort(key=lambda x: (x["chrom"], x["strand"], x["family"], x["start"], x["end"]))
with Path(snakemake.output.hsps).open("w") as out:
    for h in hsps:
        name = f'{h["family"]}|{h["query"]}|E={h["evalue"]}|bits={h["bitscore"]}'
        out.write(f'{h["chrom"]}\t{h["start"]}\t{h["end"]}\t{name}\t.\t{h["strand"]}\n')

groups = defaultdict(list)
for h in hsps:
    groups[(h["chrom"], h["strand"], h["family"])].append(h)

blast_clusters = []
for (chrom, strand, family), records in groups.items():
    current = None
    for h in records:
        if current is None or h["start"] - current["end"] > cluster_distance:
            if current is not None:
                blast_clusters.append(current)
            current = {"chrom": chrom, "start": h["start"], "end": h["end"],
                       "strand": strand, "family": family, "queries": {h["query"]}, "hsp_count": 1}
        else:
            current["end"] = max(current["end"], h["end"])
            current["queries"].add(h["query"])
            current["hsp_count"] += 1
    if current is not None:
        blast_clusters.append(current)

models = []
with miniprot_path.open() as handle:
    for line in handle:
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        name_parts = f[3].split("|", 2)
        model_id = name_parts[0]
        family = name_parts[1] if len(name_parts) > 1 else "unknown"
        models.append({"chrom": f[0], "start": int(f[1]), "end": int(f[2]),
                       "name": model_id, "family": family, "strand": f[5]})

def same_family(a, b):
    return a == b or a == "unknown" or b == "unknown"

used_models = set()
candidates = []
for cluster in blast_clusters:
    left = max(0, cluster["start"] - overlap_flank)
    right = min(chrom_lengths.get(cluster["chrom"], cluster["end"] + overlap_flank),
                cluster["end"] + overlap_flank)
    overlaps = []
    for idx, model in enumerate(models):
        if model["chrom"] != cluster["chrom"] or model["strand"] != cluster["strand"]:
            continue
        if not same_family(cluster["family"], model["family"]):
            continue
        if model["end"] > left and model["start"] < right:
            overlaps.append((idx, model))
            used_models.add(idx)
    start = min([cluster["start"]] + [m["start"] for _, m in overlaps])
    end = max([cluster["end"]] + [m["end"] for _, m in overlaps])
    candidates.append({
        "chrom": cluster["chrom"], "start": start, "end": end,
        "strand": cluster["strand"], "family": cluster["family"],
        "blast": True, "hsp_count": cluster["hsp_count"],
        "queries": sorted(cluster["queries"]), "models": [m["name"] for _, m in overlaps],
        "class": "BLAST_AND_MINIPROT" if overlaps else "BLAST_ONLY"
    })

for idx, model in enumerate(models):
    if idx in used_models:
        continue
    candidates.append({
        "chrom": model["chrom"], "start": model["start"], "end": model["end"],
        "strand": model["strand"], "family": model["family"],
        "blast": False, "hsp_count": 0, "queries": [], "models": [model["name"]],
        "class": "MINIPROT_ONLY"
    })

candidates.sort(key=lambda x: (x["chrom"], x["start"], x["end"], x["family"]))

headers = ["candidate_id", "chromosome", "start_1based", "end_1based", "strand", "family",
           "tblastn_support", "tblastn_hsp_count", "tblastn_queries", "miniprot_support",
           "miniprot_model_ids", "evidence_class"]

bed_handles = {
    "all": Path(snakemake.output.candidates).open("w"),
    "BLAST_ONLY": Path(snakemake.output.blast_only).open("w"),
    "MINIPROT_ONLY": Path(snakemake.output.miniprot_only).open("w"),
    "BLAST_AND_MINIPROT": Path(snakemake.output.both).open("w"),
}

with Path(snakemake.output.table).open("w") as table:
    table.write("\t".join(headers) + "\n")
    for number, c in enumerate(candidates, 1):
        cid = f"MHC_candidate_{number:04d}"
        row = [cid, c["chrom"], str(c["start"] + 1), str(c["end"]), c["strand"], c["family"],
               "yes" if c["blast"] else "no", str(c["hsp_count"]), ",".join(c["queries"]) or ".",
               "yes" if c["models"] else "no", ",".join(c["models"]) or ".", c["class"]]
        table.write("\t".join(row) + "\n")
        bed_line = f'{c["chrom"]}\t{c["start"]}\t{c["end"]}\t{cid}|{c["family"]}|{c["class"]}\t.\t{c["strand"]}\n'
        bed_handles["all"].write(bed_line)
        bed_handles[c["class"]].write(bed_line)

for handle in bed_handles.values():
    handle.close()

