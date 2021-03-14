nextflow.enable.dsl = 2

params.pattern = "$baseDir/data/fastq/**/*.gz"
params.executor = "local"
params.fastqc_version = "0.11.9"
params.trimmomatic_version = "0.39"
params.trimmomatic_cpus = 6

process fastqc {
  container "bschiffthaler/fastqc:${params.fastqc_version}"
  publishDir "report/00-qc/${stage}"
  executor = params.executor
  cpus = 2

  input:
    val(stage)
    tuple val(name), file(reads)

  output:
    file "*.{zip,html}"

  """
  fastqc -t 2 -o . ${reads}
  """
}

process trimmomatic {
  container "bschiffthaler/trimmomatic:${params.trimmomatic_version}"
  publishDir "analysis/01-trimmomatic"
  cpus = params.trimmomatic_cpus

  input:
    tuple val(name), file(reads)

  output:
    file "*.fastq.gz"

  """
  run.sh PE -threads ${params.trimmomatic_cpus} \
    ${reads[0]} ${reads[1]} \
    ${name}_trimmed_1.fastq.gz ${name}_orphans_1.fastq.gz \
    ${name}_trimmed_2.fastq.gz ${name}_orphans_2.fastq.gz \
    HEADCROP:5 \
    MINLEN:100
  """
}

workflow {
  data = Channel
  .fromFilePairs("data/**/P*R{1,2}*.f{astq,q}.gz", checkIfExists: true)
  
  fastqc("raw", data)
  trimmomatic(data)
}