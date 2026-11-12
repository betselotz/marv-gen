#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
================================================================================
 MARV-GEN v1.0.0 — Marburg Virus Illumina Pipeline
================================================================================
Usage:
  nextflow run main.nf -profile conda \\
    --input <fastq_dir_or_samplesheet> \\
    --host_index <bowtie2_prefix> \\
    --outdir results

Key options:
  --input              FASTQ directory or samplesheet CSV
  --host_index         Bowtie2 host index PREFIX (no .bt2)
  --reference          Marburg reference FASTA
  --reference_gff      GFF3 annotation
  --metadata           CSV with sample,date columns (TreeTime/BEAST)
  --kraken2_db         Path to existing Kraken2 DB
  --beast_xml          BEAUti XML (required for BEAST MCMC)
  --skip_beast         Skip BEAST (default: true)
  --skip_kraken2       Skip Kraken2
  --skip_host_removal  Skip host removal
  --skip_nextclade     Skip Nextclade
  --skip_phylogeny     Skip MSA/IQ-TREE
  --skip_treetime      Skip TreeTime
  --skip_gp_analysis   Skip GP analysis
================================================================================
""".stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

def input_path = params.input ?: params.samplesheet
if (!input_path) {
    log.error "Provide --input (FASTQ dir or samplesheet)"
    exit 1
}

include { INPUT_CHECK          } from './modules/input_check'
include { FASTP                } from './modules/fastp'
include { HOST_REMOVAL         } from './modules/host_removal'
include { KRAKEN2              } from './modules/kraken2'
include { KRAKEN2_BUILD        } from './modules/kraken2_build'
include { MAPPING              } from './modules/mapping'
include { PRIMER_TRIM          } from './modules/primer_trim'
include { COVERAGE_GATE        } from './modules/coverage_gate'
include { QUALIMAP             } from './modules/qualimap'
include { VARIANT_CALLING      } from './modules/variant_calling'
include { GP_ANALYSIS; GP_TREE } from './modules/gp_analysis'
include { PROTEIN_ANNOTATION   } from './modules/protein_annotation'
include { CONSENSUS            } from './modules/consensus'
include { COVERAGE             } from './modules/coverage'
include { SAMPLE_QC            } from './modules/sample_qc'
include { ANNOTATE_VARIANTS    } from './modules/annotate_variants'
include { INTRA_HOST_DIVERSITY } from './modules/intra_host_diversity'
include { NEXTCLADE_DATASET    } from './modules/nextclade_dataset'
include { NEXTCLADE            } from './modules/nextclade'
include { LINEAGE_SNPS         } from './modules/lineage_snps'
include { SNP_DISTS            } from './modules/snp_dists'
include { DN_DS                } from './modules/dn_ds'
include { COVERAGE_PLOTS       } from './modules/coverage_plots'
include { MSA                  } from './modules/msa'
include { IQTREE               } from './modules/iqtree'
include { TREETIME             } from './modules/treetime'
include { PREPARE_BEAST; BEAST_RUN; BEAST_POSTPROCESS } from './modules/beast'
include { FOCUSED_TREETIME } from './modules/focused_treetime'
include { MULTIQC              } from './modules/multiqc'

workflow {
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    ch_reference = Channel.fromPath(params.reference, checkIfExists: true).collect()
    ch_gff = file(params.reference_gff).exists() ?
        Channel.fromPath(params.reference_gff).collect() :
        Channel.empty()

    // ---- INPUT ----
    INPUT_CHECK(file(input_path))
    ch_reads = INPUT_CHECK.out.reads

    // ---- 1. FASTP ----
    FASTP(ch_reads)
    ch_trimmed = FASTP.out.reads
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(FASTP.out.json)

    // ---- 2. HOST REMOVAL ----
    if (!params.skip_host_removal) {
        ch_host = Channel.value(params.host_index.toString())
        HOST_REMOVAL(ch_trimmed, ch_host)
        ch_clean = HOST_REMOVAL.out.reads.map { meta, r1, r2 -> [meta, [r1, r2]] }
        ch_versions = ch_versions.mix(HOST_REMOVAL.out.versions.first())
    } else {
        ch_clean = ch_trimmed
    }

    // ---- KRAKEN2 ----
    if (!params.skip_kraken2) {
        def default_db = file("${projectDir}/database/kraken2_viral")
        def user_db = params.kraken2_db ? file(params.kraken2_db) : null
        if (user_db && file("${user_db}/hash.k2d").exists() && !params.build_kraken2_db) {
            ch_kdb = Channel.value(user_db)
            log.info "Kraken2 DB (user): ${user_db}"
        } else if (file("${default_db}/hash.k2d").exists() && !params.build_kraken2_db) {
            ch_kdb = Channel.value(default_db)
            log.info "Kraken2 DB (local): ${default_db}"
        } else {
            log.info "Kraken2 DB not found — building/downloading into database/kraken2_viral"
            KRAKEN2_BUILD()
            ch_kdb = KRAKEN2_BUILD.out.db
            ch_versions = ch_versions.mix(KRAKEN2_BUILD.out.versions)
        }
        KRAKEN2(ch_clean, ch_kdb)
        ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())
    }

    // ---- 3. MAPPING ----
    MAPPING(ch_clean, ch_reference)
    ch_bam = MAPPING.out.bam
    ch_versions = ch_versions.mix(MAPPING.out.versions.first())

    // ---- 3b. PRIMER TRIM (optional) ----
    if (!params.skip_primer_trim && params.primer_bed) {
        PRIMER_TRIM(ch_bam, file(params.primer_bed))
        ch_bam = PRIMER_TRIM.out.bam
        ch_versions = ch_versions.mix(PRIMER_TRIM.out.versions.first())
    }

    // ---- 4. QUALIMAP ----
    if (!params.skip_qualimap) {
        QUALIMAP(ch_bam)
        ch_versions = ch_versions.mix(QUALIMAP.out.versions.first())
    }

    // ---- 5. COVERAGE GATE ----
    COVERAGE_GATE(ch_bam)
    ch_versions = ch_versions.mix(COVERAGE_GATE.out.versions.first())
    // emit:all = [meta, bam, bai, gate.tsv] — keep PASS only as [meta, bam, bai]
    ch_bam_pass = COVERAGE_GATE.out.all
        .map { meta, bam, bai, tsv ->
            def status = "PASS"
            def f = file(tsv)
            if (f.exists()) {
                def lines = f.text.readLines()
                if (lines.size() > 1) {
                    def cols = lines[1].split('\\t')
                    status = cols[-1].trim()
                }
            }
            [meta, bam, bai, status]
        }
        .filter { meta, bam, bai, status -> status == "PASS" }
        .map { meta, bam, bai, status -> [meta, bam, bai] }

    // ---- 6. VARIANT CALLING ----
    VARIANT_CALLING(ch_bam_pass, ch_reference)
    ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions.first())

    // ---- 6b. ANNOTATE ----
    if (!params.skip_annotation) {
        ANNOTATE_VARIANTS(
            VARIANT_CALLING.out.vcf,
            file(params.reference_gff),
            file(params.reference)
        )
        ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions.first())
    }

    // ---- 6c. INTRA-HOST ----
    if (!params.skip_diversity) {
        INTRA_HOST_DIVERSITY(VARIANT_CALLING.out.tsv)
    // ch_versions = ch_versions.mix(INTRA_HOST_DIVERSITY.out.versions.first())  // process has no versions emit
    }

    // ---- 7. CONSENSUS ----
    CONSENSUS(ch_bam_pass, ch_reference)
    ch_consensus = CONSENSUS.out.fasta
    ch_versions = ch_versions.mix(CONSENSUS.out.versions.first())

    // ---- 7b. PROTEIN ANNOTATION ----
    if (!params.skip_protein_annotation && file(params.reference_gff).exists()) {
        PROTEIN_ANNOTATION(ch_consensus, file(params.reference_gff), file(params.reference))
        ch_versions = ch_versions.mix(PROTEIN_ANNOTATION.out.versions.first())

        if (!params.skip_gp_analysis) {
            def gp_comp = (params.gp_comparators_faa && file(params.gp_comparators_faa).exists()) ?
                file(params.gp_comparators_faa) : file("${projectDir}/assets/NO_FILE")
            def gp_dom = (params.gp_domains_csv && file(params.gp_domains_csv).exists()) ?
                file(params.gp_domains_csv) : file("${projectDir}/assets/GP_domains_Musoke.csv")
            GP_ANALYSIS(
                PROTEIN_ANNOTATION.out.proteins.map{ m,f -> f }.collect(),
                gp_dom,
                file(params.reference),
                file(params.reference_gff),
                gp_comp
            )
            ch_versions = ch_versions.mix(GP_ANALYSIS.out.versions)
            if (!params.skip_gp_tree) {
                GP_TREE(GP_ANALYSIS.out.alignment_faa)
                ch_versions = ch_versions.mix(GP_TREE.out.versions)
            }
        }
    }

    // ---- 8. COVERAGE ----
    ch_cov_in = ch_bam_pass.join(ch_consensus)
    COVERAGE(ch_cov_in)
    ch_versions = ch_versions.mix(COVERAGE.out.versions.first())

    if (!params.skip_sample_qc) {
        SAMPLE_QC(COVERAGE.out.stats.collect())
    }
    if (!params.skip_coverage_plots) {
        COVERAGE_PLOTS(COVERAGE.out.depth.collect())
        ch_versions = ch_versions.mix(COVERAGE_PLOTS.out.versions)
    }

    // ---- 9. NEXTCLADE ----
    if (!params.skip_nextclade) {
        NEXTCLADE_DATASET()
        ch_versions = ch_versions.mix(NEXTCLADE_DATASET.out.versions)
        NEXTCLADE(ch_consensus, NEXTCLADE_DATASET.out.dataset)
        ch_versions = ch_versions.mix(NEXTCLADE.out.versions.first())
    }

    // ---- 10. LINEAGE / dN/dS ----
    if (!params.skip_lineage_snps) {
        LINEAGE_SNPS(
            ch_consensus.map{ m,f -> f }.collect(),
            file(params.reference),
            file(params.reference_gff)
        )
        ch_versions = ch_versions.mix(LINEAGE_SNPS.out.versions)
    }
    if (!params.skip_dn_ds) {
        DN_DS(
            ch_consensus.map{ m,f -> f }.collect(),
            file(params.reference),
            file(params.reference_gff)
        )
        ch_versions = ch_versions.mix(DN_DS.out.versions)
    }

    // ---- 11–12. PHYLOGENY ----
    if (!params.skip_phylogeny) {
        ch_compare = file(params.compare_fasta_dir).exists() ?
            Channel.fromPath("${params.compare_fasta_dir}/*.fasta").collect() :
            Channel.fromPath("${projectDir}/assets/marburg_ref.fasta").collect()

        MSA(
            ch_consensus.map{ m,f -> f }.collect(),
            file(params.reference),
            ch_compare
        )
        ch_versions = ch_versions.mix(MSA.out.versions)

        IQTREE(MSA.out.alignment)
        ch_versions = ch_versions.mix(IQTREE.out.versions)

        if (!params.skip_snp_dists) {
            SNP_DISTS(MSA.out.alignment)
            ch_versions = ch_versions.mix(SNP_DISTS.out.versions)
        }

        if (!params.skip_treetime && params.metadata) {
            TREETIME(
                IQTREE.out.tree,
                MSA.out.alignment,
                file(params.metadata)
            )
            ch_versions = ch_versions.mix(TREETIME.out.versions)
        }

        // Focused TreeTime (JN408064.1 + Ethiopian sequences only)
        if (params.focused_treetime && params.metadata) {
            FOCUSED_TREETIME(
                MSA.out.alignment,
                file(params.metadata)
            )
            ch_versions = ch_versions.mix(FOCUSED_TREETIME.out.versions)
        }

        // BEAST: dates always if enabled; MCMC only with --beast_xml
        if (!params.skip_beast) {
            def meta_file = params.metadata ? file(params.metadata) : file("${projectDir}/metadata/all_seq_metadata.csv")
            if (meta_file.exists()) {
                PREPARE_BEAST(MSA.out.alignment, meta_file)
                ch_versions = ch_versions.mix(PREPARE_BEAST.out.versions)

                def xml_path = params.beast_xml ? file(params.beast_xml.toString()) : null
                if (xml_path && xml_path.exists()) {
                    ch_xml = Channel.value(file(params.beast_xml))
                    ch_run = Channel.of(
                        [1, (params.beast_seed1 ?: 1111)],
                        [2, (params.beast_seed2 ?: 2222)]
                    )
                    BEAST_RUN(ch_xml, ch_run.map{ it[0] }, ch_run.map{ it[1] })
                    ch_logs  = BEAST_RUN.out.log.toSortedList{ a, b -> a[0] <=> b[0] }
                    ch_trees = BEAST_RUN.out.trees.toSortedList{ a, b -> a[0] <=> b[0] }
                    BEAST_POSTPROCESS(
                        ch_logs.map{ it[0][1] },
                        ch_logs.map{ it[1][1] },
                        ch_trees.map{ it[0][1] },
                        ch_trees.map{ it[1][1] }
                    )
                    ch_versions = ch_versions.mix(BEAST_RUN.out.versions)
                    ch_versions = ch_versions.mix(BEAST_POSTPROCESS.out.versions)
                } else {
                    log.warn "BEAST: no --beast_xml — only tip dates/alignment prepared (see scripts/beast/02_generate_beast_xml.sh)"
                }
            } else {
                log.warn "BEAST enabled but metadata not found"
            }
        }
    }

    // ---- MULTIQC ----
    if (!params.skip_multiqc) {
        MULTIQC(Channel.value(file(params.outdir)))
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    ch_versions.unique().collectFile(
        name: 'software_versions.yml',
        storeDir: "${params.outdir}/pipeline_info"
    )
}

workflow.onComplete {
    log.info """
================================================================================
 MARV-GEN v1.0.0 completed
 Status  : ${workflow.success ? 'SUCCESS' : 'FAILED'}
 Duration: ${workflow.duration}
 Outdir  : ${params.outdir}
================================================================================
""".stripIndent()
}
