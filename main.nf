
log.info """\
    RNA-SEQ  RIBO-DEPLETED  PIPELINE
    ================================
    genome       : ${params.genomeFasta}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()


// ====================
// FastQC process
// ====================

process FastQC {

  publishDir "${params.outdir}/qc", mode: 'copy'

  input:
    path fastq_file

  output:
    path "*.html", emit: qc_reports
    path "*.zip",  emit: qc_zips


  script:
  
  sample_id = fastq_file.getBaseName().replaceAll(/\.fastq\.gz$/, '')
  """
  fastqc ${fastq_file}
  """
}

// ====================
// MultiQC process
// ====================

process MultiQC {

  tag "MultiQC summary"

  input:
  path fastqc_results

  output:
  path "multiqc_report.html",  emit: multiqc_html
  path "multiqc_data",         emit: multiqc_data

  publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

  script:
  """
  multiqc ${fastqc_results} --outdir . 
  """
}

// ====================
// Trimmomatic process
// ====================

process Trimmomatic {
  
  publishDir "${params.outdir}/trimmed", mode: "copy"

  input:
  tuple val(sample_id), path(reads)

  output: 
  tuple val(sample_id), 
        path("${sample_id}_R1_trimmed.fastq.gz"),
        path("${sample_id}_R2_trimmed.fastq.gz"),
        emit: trimmed_pairs
  
  script:
  adapters = "${CONDA_PREFIX}/share/trimmomatic-*/adapters/TruSeq3-PE.fa"
  """
  trimmomatic PE -threads ${task.cpus} -phred33 \
    ${reads[0]} ${reads[1]} \
    ${sample_id}_R1_trimmed.fastq.gz /dev/null \
    ${sample_id}_R2_trimmed.fastq.gz /dev/null \
    ILLUMINACLIP:${adapters}:2:30:10:2:keepBothReads \
    SLIDINGWINDOW:5:20 \
    MINLEN:17

  """
}

// ====================
// STAR process
// ====================

process STARindex {
  tag "STAR index"
  publishDir params.starIndex, mode: 'copy', overwrite: false

  input:
  path genomeFasta
  path gtf

  output:
  path "star_index", emit: star_index_dir

  script:
  """
  # Decompress inputs
  gunzip -c ${genomeFasta} > genome.fa
  gunzip -c ${gtf}         > annotation.gtf

  # Build the STAR index
  mkdir -p star_index
  STAR \
    --runThreadN ${task.cpus} \
    --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles genome.fa \
    --sjdbGTFfile annotation.gtf \
    --sjdbOverhang 58

  """
}

process STAR {
  tag "${sample_id}"

  // Publish BAMs to alignment/
  publishDir "${params.outdir}/alignment",
             mode: 'copy',
             pattern: "*Aligned.sortedByCoord.out.bam"

  // Publish gene counts to counts/
  publishDir "${params.outdir}/counts",
             mode: 'copy',
             pattern: "*ReadsPerGene.out.tab"

  input:
    tuple val(sample_id), path(r1), path(r2)
    path star_index
    path gtf

  output:
    tuple val(sample_id),
          path("${sample_id}Aligned.sortedByCoord.out.bam"),
          emit: aligned_bam
    tuple val(sample_id),
          path("${sample_id}ReadsPerGene.out.tab"),
          emit: gene_counts 


  script:
    """
    gunzip -c ${gtf} > annotation.gtf

    STAR \
      --quantMode GeneCounts \
      --runThreadN ${task.cpus} \
      --genomeDir ${star_index} \
      --sjdbGTFfile annotation.gtf \
      --readFilesIn ${r1} ${r2} \
      --readFilesCommand zcat \
      --outFileNamePrefix ${sample_id} \
      --outSAMtype BAM SortedByCoordinate
    """
}


// ====================
// Workflow
// ====================

workflow {

  /* 1. raw paired-end reads */
  raw_reads = Channel
                .fromFilePairs( params.reads )         // emits [ id , [R1,R2] ]

  /* 2. trim */
  Trimmomatic( raw_reads )

  /* 3. QC on trimmed reads */
  trimmed_fastq = Trimmomatic.out.trimmed_pairs
                    .map { id, r1, r2 -> [ r1, r2 ] }  // drop id, keep paths
                    .flatten()                         // one path per item
  raw_fastq = raw_reads
                  .map{ id, pair -> pair }             // give FastQC the paths only
                  .flatten()
  FastQC( trimmed_fastq.mix(raw_fastq) )

  /* 5. collect FastQC ZIPs (raw + trimmed) for MultiQC */
  qc_zips = FastQC.out.qc_zips.collect()
  MultiQC( qc_zips )

  /* 6. Build STAR index */
  genome_fa  = Channel.value( file(params.genomeFasta) )
  genome_gtf = Channel.value( file(params.gtf) )
  star_idx   = STARindex(genome_fa, genome_gtf).star_index_dir

  /* 7. Align trimmed reads */
  alignment = STAR(Trimmomatic.out.trimmed_pairs, 
                   star_idx,
                   genome_gtf)
}
