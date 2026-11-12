# MARV-GEN

**MARV-GEN** is a Nextflow DSL2 pipeline for Illumina whole-genome sequencing (WGS) analysis of Marburg virus (MARV).

Developed by the Ethiopian Public Health Institute (EPHI) Bioinformatics team.

The pipeline integrates read quality control, host-read removal, reference mapping, variant calling, consensus generation, genomic characterization, lineage assignment, and phylogenetic analysis.

## Features

- Read quality control and trimming with fastp
- Host-read removal using a Bowtie2 index
- Optional taxonomic classification using Kraken2
- Reference mapping and variant calling
- Consensus genome generation
- Coverage assessment and sample-level quality control
- Nextclade lineage and clade assignment
- Variant annotation and protein-level analysis
- Glycoprotein (GP) analysis and comparison with reference sequences
- Multiple sequence alignment and phylogenetic inference
- SNP distance analysis
- Temporal phylogenetic analysis using TreeTime
- Bayesian phylogenetic analysis using BEAST
- MultiQC reporting

## Requirements

| Requirement | Version / information |
|---|---|
| Operating system | Linux |
| Nextflow | >=23.04.0 |
| Java | OpenJDK 17 recommended |
| Conda | Miniconda or Anaconda |
| Memory | 16 GB or more recommended |
| CPUs | 8 or more recommended |

The pipeline uses Conda environments defined in `envs/` to manage process-level software dependencies.

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/betselotz/marv-gen.git
cd marv-gen
```

### 2. Install dependencies

```bash
bash bin/marv-gen-setup
```

The setup script creates a Conda environment named `nf-java17`, installs Java 17 and Nextflow, installs the MARV-GEN command-line wrapper, and prepares the process environments.

Environment creation may take some time.

### 3. Activate the environment

```bash
source ~/.bashrc
conda activate nf-java17
```

Check the installation:

```bash
marv-gen --help
nextflow -version
java -version
```

## Input preparation

MARV-GEN accepts either:

1. A CSV samplesheet containing paired-end FASTQ file paths.
2. A directory containing paired-end FASTQ files named using the `_R1` and `_R2` convention.

### Samplesheet format

Prepare a CSV file with three columns:

| sample | fastq_1 | fastq_2 |
|---|---|---|
| MARV_001 | /data/MARV_001_R1.fastq.gz | /data/MARV_001_R2.fastq.gz |
| MARV_002 | /data/MARV_002_R1.fastq.gz | /data/MARV_002_R2.fastq.gz |

An example is provided in:

`assets/samplesheet.example.csv`

Use absolute paths to FASTQ files. Sample identifiers should be unique.

## Reference preparation

### Marburg reference genome

The default reference genome and annotation are included:

- `reference_genomes/Marburg_reference.fasta`
- `reference_genomes/Marburg_reference.gff3`

Comparison sequences are available in:

`reference_genomes/MARV_compare/`

### Host-read removal

Host-read removal is enabled by default and requires a Bowtie2 index of the human reference genome.

Obtain a human reference genome FASTA, for example GRCh38, and build the index:

```bash
conda install -c bioconda bowtie2 -y

bowtie2-build GRCh38.fa reference_genomes/human_index
```

The resulting index prefix is:

`reference_genomes/human_index`

Alternatively, supply another index prefix using `--host_index`.

If host-read removal is not required, it can be skipped:

```bash
marv-gen --input samplesheet.csv --skip_host_removal true
```

## Running the pipeline

### Basic analysis

With a prepared human Bowtie2 index:

```bash
marv-gen --input samplesheet.csv
```

To run using a directory of paired-end FASTQ files:

```bash
marv-gen --input /data/marv_fastq/
```

To resume an interrupted analysis:

```bash
marv-gen --input samplesheet.csv -resume
```

To specify an output directory:

```bash
marv-gen --input samplesheet.csv --outdir results
```

### Adjusting quality-control parameters

Default consensus quality thresholds include:

- Minimum coverage percentage: 70
- Maximum ambiguous nucleotide percentage (N): 20

Example:

```bash
marv-gen --input samplesheet.csv \
  --min_coverage_pct 50 \
  --max_n_pct 25
```

### Optional analyses

The following components are disabled by default and can be enabled through command-line parameters.

Enable TreeTime:

```bash
marv-gen --input samplesheet.csv --skip_treetime false
```

Enable BEAST:

```bash
marv-gen --input samplesheet.csv --skip_beast false
```

Enable glycoprotein analysis:

```bash
marv-gen --input samplesheet.csv --skip_gp_analysis false
```

Enable Kraken2:

```bash
marv-gen --input samplesheet.csv \
  --skip_kraken2 false \
  --kraken2_db /path/to/kraken2_database
```

BEAST analysis uses the XML configuration provided in `assets/`. A custom XML file can be supplied using `--beast_xml`.

## Default configuration

| Component | Default |
|---|---|
| Host-read removal | Enabled |
| Kraken2 | Disabled |
| Primer trimming | Disabled |
| Nextclade | Enabled |
| Phylogenetic analysis | Enabled |
| TreeTime | Disabled |
| BEAST | Disabled |
| GP analysis | Disabled |
| MultiQC | Enabled |

These settings can be changed using command-line parameters.

## Workflow outputs

Results are written to the directory specified by `--outdir` (default: `results/`).

| Directory | Description |
|---|---|
| `02_nonhuman_reads/` | Reads retained after host-read removal |
| `07_consensus/` | Consensus genome sequences |
| `08_coverage/` | Coverage summaries and sample QC |
| `09_nextclade/` | Nextclade lineage and clade assignments |
| `10_lineage/` | Lineage-associated SNP and related analyses |
| `11_phylogeny/` | Multiple sequence alignment, IQ-TREE and SNP distance outputs |
| `12_focused_treetime/` | Focused temporal phylogenetic analysis |
| `12_timetree/` | Temporal phylogenetic outputs, including BEAST when enabled |
| `13_multiqc/` | MultiQC reports |
| `pipeline_info/` | Pipeline execution and software-version information |

Additional outputs may be generated by optional analysis modules.

## Metadata and temporal analysis

Sequence metadata are provided in:

`metadata/all_seq_metadata.csv`

For TreeTime and BEAST analyses, ensure that sequence identifiers in the metadata match the identifiers in the alignment and that collection dates are correctly formatted.

## Reproducibility

For reproducible analyses:

- Record the MARV-GEN release and Nextflow version.
- Retain the original samplesheet and input sequence information.
- Record the reference genome and database versions.
- Preserve the software-version report generated in `pipeline_info/`.
- Report the parameters and optional modules used in the analysis.
- Use `-resume` when restarting an interrupted Nextflow execution.

Default quality-control thresholds and analysis settings are defined in `nextflow.config` and can be overridden from the command line.

## Repository structure

```text
marv-gen/
├── assets/                 # Example inputs and analysis configuration
├── bin/                    # Command-line and setup scripts
├── conf/                   # Pipeline configuration
├── database/               # MARV Nextclade dataset
├── envs/                   # Conda environment definitions
├── metadata/               # Sequence metadata
├── modules/                # Nextflow process modules
├── reference_genomes/      # MARV reference and comparison sequences
├── main.nf                 # Main Nextflow workflow
├── nextflow.config         # Pipeline configuration
├── install.sh
└── README.md
```

## Citation

If you use MARV-GEN in your research, please cite the repository and the specific release used:

**MARV-GEN: Marburg Virus Illumina Genome Analysis Pipeline.** EPHI Bioinformatics.

https://github.com/betselotz/marv-gen

A peer-reviewed publication describing the pipeline will be added when available.

## Acknowledgements

Developed by the Ethiopian Public Health Institute (EPHI) Bioinformatics team, supporting pathogen genomics and public health surveillance.

## License

See the `LICENSE` file for licensing information.
