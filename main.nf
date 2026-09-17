#!/usr/bin/env nextflow
/*
 * CRAM -> BAM conversion pipeline.
 *
 * CRAM files reference-compress alignment data and cannot be decoded
 * without the exact reference FASTA (+ .fai) used at compression time --
 * unlike BAM, which is self-contained. samtools view handles the
 * decompression and format conversion in one step.
 *
 * Input: params.cram_runs, injected by preprocess.py from the Cirro
 * dataset's file listing -- NOT a user-supplied samplesheet.
 * Output: <outdir>/<sample_id>/<sample_id>.bam (+ .bam.bai)
 *         <outdir>/samplesheet.csv -- columns: sample,file
 *         (one row per BAM, one row per BAI, so each sample_id appears twice)
 */
nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// PARAMS
// ---------------------------------------------------------------------------

if (!params.containsKey('outdir'))     params.outdir = 'results'
if (!params.containsKey('ref_fasta'))  params.ref_fasta = null
if (!params.containsKey('ref_fai'))    params.ref_fai = null

// ---------------------------------------------------------------------------
// PROCESSES
// ---------------------------------------------------------------------------

process CRAM_TO_BAM {
    tag "${sample_id}"
    label 'process_medium'
    container "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
    errorStrategy 'retry'
    maxRetries 3

    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(cram), path(crai)
    path ref_fasta
    path ref_fai

    output:
    tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), emit: bam

    shell:
    '''
    set -euxo pipefail

    if [ ! -f !{cram}.crai ]; then
        ln -s !{crai} !{cram}.crai
    fi

    samtools view \
        -@ !{task.cpus} \
        -b \
        -T !{ref_fasta} \
        -o !{sample_id}.bam \
        !{cram}

    samtools index -@ !{task.cpus} !{sample_id}.bam

    test -s !{sample_id}.bam
    test -s !{sample_id}.bam.bai
    '''
}

process WRITE_SAMPLESHEET {
    tag "samplesheet"
    label 'process_low'
    container "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
    publishDir "${params.outdir}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2

    input:
    val csv_content

    output:
    path "samplesheet.csv", emit: samplesheet

    shell:
    '''
    cat > samplesheet.csv << 'CSV_EOF'
!{csv_content}
CSV_EOF

    test -s samplesheet.csv
    '''
}

// ---------------------------------------------------------------------------
// WORKFLOW
// ---------------------------------------------------------------------------

workflow {

    if (!params.cram_runs) {
        error "params.cram_runs is empty -- did the Cirro preprocess.py hook run? (see preprocess.py)"
    }
    if (!params.ref_fasta) {
        error "Missing required param: ref_fasta"
    }
    if (!params.ref_fai) {
        error "Missing required param: ref_fai"
    }

    ref_fasta = file(params.ref_fasta, checkIfExists: true)
    ref_fai   = file(params.ref_fai,   checkIfExists: true)

    samples_ch =
        Channel
            .fromList(params.cram_runs)
            .map { run ->
                tuple(
                    run.sample_id,
                    file(run.cram, checkIfExists: true),
                    file(run.crai, checkIfExists: true)
                )
            }

    CRAM_TO_BAM(
        samples_ch,
        ref_fasta,
        ref_fai
    )

    // -----------------------------------------------------------------------
    // Collect every sample's [sample_id, bam, bai] into one list, then emit
    // the samplesheet once all conversions have finished -- .collect()
    // forces this to wait for every CRAM_TO_BAM call to complete first.
    // -----------------------------------------------------------------------

    samplesheet_rows =
    CRAM_TO_BAM
        .out
        .bam
        .map { sample_id, bam, bai -> [sample_id, bam, bai] }
        .collect(flatten: false)

    samplesheet_content =
        samplesheet_rows.map { rows ->
            def lines = ["sample,file"]
            rows.each { sample_id, bam, bai ->
                lines << "${sample_id},${params.outdir}/${sample_id}/${sample_id}.bam"
                lines << "${sample_id},${params.outdir}/${sample_id}/${sample_id}.bam.bai"
            }
            lines.join("\n")
        }

    WRITE_SAMPLESHEET(samplesheet_content)
}