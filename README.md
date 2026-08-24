# Avian MHC reference-annotation pipeline

A reproducible discovery and evidence-integration workflow for annotating major
histocompatibility complex (MHC) genes in bird reference genomes.

The workflow supports genomes **with or without an existing annotation**:

1. **BLASTP** audits MHC proteins already present in an annotated proteome.
2. **TBLASTN** searches the genomic assembly independently for MHC-like coding sequence.
3. **miniprot** constructs splice-aware protein-to-genome models.
4. A reconciliation script groups TBLASTN HSPs into candidate loci and compares them
   with miniprot models.
5. The resulting evidence package is reviewed manually and converted into a curated GFF3.

> This pipeline discovers and organizes evidence. It deliberately does not assign
> definitive locus names or automatically declare pseudogenes. Avian MHC regions are
> duplicated, polymorphic, and assembly-sensitive; final annotation requires biological
> review, RNA/genomic-read evidence, phylogeny, and synteny.

## Workflow

```text
Related-bird MHC proteins
       |                         Existing annotation (optional)
       |                                      |
       |                                  BLASTP audit
       |                                      |
       +---------- TBLASTN genome search -----+
       |                                      |
       +---------- miniprot gene models ------+
                              |
                coordinate reconciliation
                              |
        candidate table + BED + evidence files
                              |
                manual validation and GFF3
```

## Why use all three searches?

| Method | Query -> target | Purpose | Important limitation |
|---|---|---|---|
| BLASTP | protein -> predicted proteins | Find MHC models already present in an annotation | Cannot find a gene omitted or mistranslated by that annotation |
| TBLASTN | protein -> genomic DNA translated in six frames | Discover coding fragments anywhere in the assembly | Reports local HSPs, often corresponding to separate exons |
| miniprot | protein -> genomic DNA | Construct splice-aware CDS models | May disfavor weak fragments, damaged pseudogenes, or very divergent copies |
| BLASTN | nucleotide -> genomic DNA | Optional search for close-species exons, CDS, alleles, or flanks | Less sensitive across evolutionary distance; not part of the default workflow |

## Repository layout

```text
avian-mhc-annotation-pipeline/
├── README.md
├── LICENSE
├── CITATION.cff
├── config/
│   └── config.example.yaml
├── envs/
│   └── mhc.yaml
├── workflow/
│   ├── Snakefile
│   └── scripts/
│       ├── miniprot_gff_to_bed.py
│       └── reconcile_candidates.py
└── examples/
    ├── reference_protein_headers.example.txt
    └── protein_metadata.example.tsv
```

## Requirements

- Linux or macOS
- Conda/Mamba
- Snakemake >= 8
- Input genome FASTA
- Full-length MHC/reference proteins in FASTA
- Optionally, an existing GFF3 and predicted proteome

The supplied Conda environment contains:

- NCBI BLAST+
- miniprot
- samtools
- gffread
- seqkit
- Python

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/avian-mhc-annotation-pipeline.git
cd avian-mhc-annotation-pipeline

mamba env create -f envs/mhc.yaml
conda activate avian-mhc-annotation
```

Alternatively, let Snakemake create rule-specific environments with `--use-conda`.

## Prepare reference proteins

Collect full-length proteins from several related bird species. At minimum, search
separate representatives of:

- MHC class I
- MHC class II alpha
- MHC class II beta
- TAP1
- TAP2
- TAPBP/tapasin
- DMA and DMB
- B2M and other pathway genes if desired

Recommended FASTA identifiers:

```text
>MHCI|Gallus_gallus|BF2|NP_xxxxx|complete
MAVMAPR...
>MHCIIB|Close_species|unknown|XP_xxxxx|predicted
MTRLL...
```

The pipeline uses the text before the first `|` as the protein family when clustering
TBLASTN hits. If headers do not contain `|`, the complete query ID becomes the family.

Keep full-length proteins separate from partial exon/allele sequences. Full-length
proteins are appropriate for miniprot. Partial sequences can be searched separately
with TBLASTN or BLASTN during manual follow-up.

## Configure a run

Copy the example configuration:

```bash
cp config/config.example.yaml config/config.yaml
```

Edit paths in `config/config.yaml`:

```yaml
genome: "/absolute/path/species.genome.fa"
reference_proteins: "/absolute/path/mhc_reference_proteins.fa"

# Optional. Set to "" when unavailable.
existing_gff: "/absolute/path/species.annotation.gff3"
annotated_proteome: "/absolute/path/species.proteins.faa"

threads: 16
blast:
  tblastn_evalue: 1e-3
  blastp_evalue: 1e-5
  max_target_seqs: 1000
miniprot:
  outn: 50
  outs: 0.5
  outc: 0.25
reconciliation:
  cluster_distance: 20000
  overlap_flank: 10000
```

Use absolute paths where possible. If `annotated_proteome` is empty but `existing_gff`
is supplied, the workflow translates proteins from the GFF3 with gffread. If both are
empty, the BLASTP branch writes an explanatory skip file and the genome-level workflow
continues normally.

## Run

Preview the directed acyclic graph:

```bash
snakemake -s workflow/Snakefile \
  --configfile config/config.yaml \
  --use-conda \
  -n
```

Execute:

```bash
snakemake -s workflow/Snakefile \
  --configfile config/config.yaml \
  --use-conda \
  --cores 16 \
  --printshellcmds
```

## Outputs

```text
results/
├── 00_qc/
│   ├── genome.fa.fai
│   └── reference_protein_lengths.tsv
├── 01_blastp/
│   ├── mhc_vs_annotated_proteome.tsv
│   └── SKIPPED.txt                  # only when no annotation/proteome is supplied
├── 02_tblastn/
│   ├── mhc_vs_genome.tsv
│   └── tblastn_hsps.bed
├── 03_miniprot/
│   ├── mhc.miniprot.gff3
│   └── miniprot_models.bed
└── 04_reconciliation/
    ├── candidate_loci.tsv
    ├── candidate_loci.bed
    ├── blast_only_candidates.bed
    ├── miniprot_only_candidates.bed
    └── supported_by_both.bed
```

### `candidate_loci.tsv`

The principal table contains one row per nonredundant candidate region:

```text
candidate_id, chromosome, start, end, strand, family,
tblastn_support, tblastn_hsp_count, tblastn_queries,
miniprot_support, miniprot_model_ids, evidence_class
```

Evidence classes are:

- `BLAST_AND_MINIPROT`
- `BLAST_ONLY`
- `MINIPROT_ONLY`

These classes describe computational evidence, not gene functionality.

## How reconciliation works

1. TBLASTN HSPs are converted from 1-based BLAST coordinates to 0-based BED.
2. HSPs are grouped by chromosome, strand, and query family.
3. Nearby HSPs within `cluster_distance` are combined into provisional candidate regions.
4. Candidate regions are expanded by `overlap_flank` for comparison with miniprot.
5. Same-strand overlaps are reported and collapsed into a nonredundant evidence table.

The clustering distance is only a discovery parameter. In a tandem MHC expansion,
two adjacent genes may be placed in one provisional region. Final gene boundaries must
come from miniprot, RNA junctions, ORFs, and manual inspection—not the BLAST cluster.

## Manual curation after the workflow

For each candidate:

1. Load TBLASTN HSPs, miniprot GFF3, the existing annotation, RNA BAMs, and genomic-read
   BAMs in IGV, JBrowse, or Apollo.
2. Confirm exon order, splice sites, CDS phase, start/stop codons, and ORF.
3. Check signal peptide, MHC domains, transmembrane helix, and cytoplasmic tail.
4. Look for unsupported fusions of adjacent paralogues.
5. Verify frameshifts and premature stops using genomic reads before calling a pseudogene.
6. Inspect contig boundaries, gaps, abnormal coverage, and alternative haplotigs.
7. Use family-specific phylogenies plus synteny before transferring locus names.

Suggested provisional names:

```text
Species_MHCI_1
Species_MHCI_2
Species_MHCIIA_1
Species_MHCIIB_1
```

## Interpreting discordant results

| Evidence | Likely interpretation | Recommended action |
|---|---|---|
| BLASTP + TBLASTN + complete miniprot | Existing annotation is probably credible | Validate ORF, domains, RNA, and synteny |
| No BLASTP; TBLASTN + complete miniprot | Existing annotation likely missed the gene | Build a new curated model |
| TBLASTN only | Partial gene, pseudogene, divergent exon, or false positive | Inspect locally; consider Exonerate/GeneWise |
| miniprot only | Weak individual exons or permissive model | Rerun local TBLASTN and check coverage/domains |
| Many HSPs but one model | Tandem duplication or collapsed copies possible | Inspect coverage and long-read evidence |
| One model spans repeating HSP patterns | Possible fusion of adjacent paralogues | Check RNA junctions and separate ORFs |

## RNA-seq

RNA-seq is not required by this discovery workflow. When available, align it separately
with STAR or HISAT2 and add splice junctions/BAM coverage during manual curation.

Useful tissues include spleen, thymus, bursa, blood/leukocytes, lung, intestine, and
infected or stimulated tissues. Absence of RNA support is not evidence that a genomic
MHC gene is absent.

## Preparing for resequencing analyses

The final manually curated annotation should provide:

- `mhc_curated.gff3`
- complete gene spans in BED
- exon and CDS BED files
- MHC-I exons 2/3 and MHC-IIB exon 2 BED files after locus-specific verification
- genes plus flanking intervals
- callable, uniquely mappable subregions
- CDS and predicted protein FASTA files
- an evidence/confidence table including assembly warnings

Map resequencing reads to the complete genome, not only extracted MHC sequences. Before
interpreting SNP diversity, inspect mapping quality, depth, secondary alignments, allele
balance, copy-number variation, and structural variation.

## Reproducibility

- Record the assembly accession and exact FASTA checksum.
- Preserve protein accessions and source species in a metadata table.
- Save the completed configuration with each run.
- Commit manual GFF3 edits and evidence decisions to version control.
- Pin software versions when producing a publication release.

## License

MIT. See [LICENSE](LICENSE).

