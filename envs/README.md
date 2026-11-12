# Conda environments used by marv-gen

| File | Used by |
|------|---------|
| focused_treetime_env.yaml | FOCUSED_TREETIME |
| python_analysis_env.yaml | PREPARE_BEAST and helpers |
| beast_env.yaml | BEAST_RUN and BEAST_POSTPROCESS |

Host env for Nextflow:

  conda create -n nf-java17 -c conda-forge -c bioconda openjdk=17 nextflow -y
  conda activate nf-java17
