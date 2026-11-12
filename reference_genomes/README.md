# Reference genomes

This directory contains reference sequences used by MARV-GEN.

## Included files

- `Marburg_reference.fasta`: MARV reference genome used for mapping and consensus generation.
- `Marburg_reference.gff3`: Reference genome annotation.
- `MARV_compare/`: MARV comparison sequences used in comparative analyses.
- `KC545387.fasta`: Additional reference sequence.

## Human reference genome for host-read removal

The human reference genome and its Bowtie2 index are not distributed with this repository.

Host-read removal is enabled by default. Obtain a suitable human reference genome FASTA and build the index:

```bash
bowtie2-build GRCh38.fa human_index
```

Supply the resulting index prefix to MARV-GEN:

```bash
marv-gen --input samplesheet.csv --host_index reference_genomes/human_index
```

Alternatively, skip host-read removal:

```bash
marv-gen --input samplesheet.csv --skip_host_removal true
```
