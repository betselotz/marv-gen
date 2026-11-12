def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id = row.sample.trim()
    def fastq_1 = file(row.fastq_1.trim())
    def fastq_2 = row.fastq_2 ? file(row.fastq_2.trim()) : null
    if (!fastq_1.exists()) exit 1, "ERROR: Read 1 missing: ${row.fastq_1}"
    if (fastq_2 && !fastq_2.exists()) exit 1, "ERROR: Read 2 missing: ${row.fastq_2}"
    return [meta, [fastq_1, fastq_2].findAll{ it }]
}

workflow INPUT_CHECK {
    take:
    input_path
    main:
    def path_str = input_path.toString()
    if (path_str.endsWith('.csv') || path_str.endsWith('.tsv')) {
        Channel
            .fromPath(input_path, checkIfExists: true)
            .splitCsv(header: true, sep: path_str.endsWith('.tsv') ? '\t' : ',')
            .map { create_fastq_channel(it) }
            .set { reads }
    } else {
        Channel
            .fromFilePairs("${path_str}/*_R{1,2}*.fastq.gz", checkIfExists: true)
            .ifEmpty { error "ERROR: No paired FASTQ files in: ${path_str}" }
            .map { id, files ->
                def sample = id.replaceAll(/_R[12]$/, '')
                [ [id: sample], files.sort() ]
            }
            .set { reads }
    }
    reads
        .map { meta, files -> meta.id }
        .collect()
        .subscribe { samples ->
            log.info "Detected ${samples.size()} sample(s): ${samples.join(', ')}"
        }
    emit:
    reads
}
