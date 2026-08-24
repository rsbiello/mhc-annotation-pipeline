from pathlib import Path

source = Path(snakemake.input[0])
target = Path(snakemake.output[0])
target.parent.mkdir(parents=True, exist_ok=True)

with source.open() as inp, target.open("w") as out:
    for line in inp:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or fields[2] != "mRNA":
            continue
        chrom, start, end, strand, attrs = fields[0], int(fields[3]), int(fields[4]), fields[6], fields[8]
        attr_map = {}
        for token in attrs.split(";"):
            if "=" in token:
                key, value = token.split("=", 1)
                attr_map[key] = value
        model_id = attr_map.get("ID", f"miniprot_{chrom}_{start}_{end}")
        query = attr_map.get("Target", attr_map.get("Query", "unknown")).split()[0]
        family = query.split("|", 1)[0] if query != "unknown" else "unknown"
        name = f"{model_id}|{family}|{query}"
        out.write(f"{chrom}\t{start-1}\t{end}\t{name}\t.\t{strand}\n")

