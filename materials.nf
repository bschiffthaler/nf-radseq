nextflow.enable.dsl=2

include { citekey as fq_citekey } from './modules/fastqc.nf'
include { bibtex as fq_bib } from './modules/fastqc.nf'

include { citekey as tr_citekey } from './modules/trimmomatic.nf'
include { bibtex as tr_bib } from './modules/trimmomatic.nf'

include { citekey as st_citekey } from './modules/stacks.nf'
include { bibtex as st_bib } from './modules/stacks.nf'

include { citekey as mq_citekey } from './modules/multiqc.nf'
include { bibtex as mq_bib } from './modules/multiqc.nf'

process gen_materials {

  container "bschiffthaler/latex"
  executor params.executor
  cpus 1
  publishDir "report/materials"

  input:
  val "dummy"

  output:
  path "materials.md"
  path "references.bib"
  path "materials.docx"

  script:

  fq_key = fq_citekey()
  fq_bib = fq_bib()

  tr_key = tr_citekey()
  tr_bib = tr_bib()
  tr_trimmers = params.trimmomatic_trimmers.join(" ")

  st_key = st_citekey()
  st_bib = st_bib()
  pr_params = params.process_radtags_params.join(" ")
  us_params = params.ustacks_params.join(" ")
  cs_params = params.cstacks_params.join(" ")
  ss_params = params.sstacks_params.join(" ")
  ts_params = params.tsv2bam_params.join(" ")
  po_params = params.populations_params.join(" ")

  mq_key = mq_citekey()
  mq_bib = mq_bib()

  """
  echo "# Materials and methods
  Unless otherwise stated all program options were left default.

  Raw and trimmed data quality were assessed using FastQC
  (${params.fastqc_version}) [${fq_key}]. The raw data was trimmed using
  trimmomatic (${params.trimmomatic_version}) in PE mode [${tr_key}].
  The following trimmers were used in order:

  ~~~~~~
  ${tr_trimmers}
  ~~~~~~

  Next, data was processed through stacks (${params.stacks_version})[${st_key}].

  ~~~~~~
  process_radtags ${pr_params}
  ustacks ${us_params}
  sstacks ${ss_params}
  cstacks ${cs_params}
  tsv2bam ${ts_params}
  populations ${po_params}
  ~~~~~~

  Finally, supported logs and data were aggregated into a report using
  MultiQC (${params.multiqc_version}) [${mq_key}].
  " > materials.md

  echo "${fq_bib}
  ${tr_bib}
  ${st_bib}
  ${mq_bib}" > references.bib

  pandoc -i materials.md -o materials.docx -C --bibliography references.bib
  """

}

workflow {
  gen_materials("")
}