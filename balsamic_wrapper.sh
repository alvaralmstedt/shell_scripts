#!/bin/bash -l
set -eo pipefail
shopt -s expand_aliases

_log="[$(date) $(whoami)] "
_red=${_log}'\033[0;31m';
_green=${_log}'\033[0;32m';
_yellow=${_log}'\033[1;33m';
_nocol='\033[0m';

function usage() {
  echo $"
USAGE: [ -a <T|TN> -c _condaenv -m <run|config|all> -t <panel|WGS> -r ]
  
  -a [required] T: Tumor only, TN: tumor normal
  -c Conda environment where BALSAMIC is installed. If not specified, it will use current environment.
  -m [required] config: only create config file, run: create config file and start analysis
  -t [required] twist_ffpe_exome Exomes for , GMS-ST: Twist GMS-ST panel, test_panel: target sequencing workflow (also includes WES), WGS: whole genome sequencing workflow
  -o [required] outout path
  -d Analysis dir path, if it doesn't exist it will be created
  -r Flag. Set to submit jobs instead of running in dry mode
  -f [required] tumor forward fastq path.Reverse should be in the same folder and wil be detected autmatically.
  -n normal fastq if TN is chosen
  -h Show this help and exit
" 
}

while getopts ":a:c:m:t:d:o:f:n:r" opt; do
  case ${opt} in
    a)
      _analysis=${OPTARG}
      echo "analysis set to" "${OPTARG}"
      [[ $_analysis == 'T' || $_analysis == 'TN' ]] || ( usage >&2; exit 1)
      ;;
    c)
      _condaenv=${OPTARG}
      echo "conda environment set to" "${OPTARG}"
      ;;
    m)
      _startmode=${OPTARG}
      echo "start mode set to" "${OPTARG}"
      [[ $_startmode == 'config' || $_startmode == 'run' || $_startmode == 'all' ]] || ( usage >&2; exit 1)
      ;;
    t)
      _ngstype=${OPTARG}
      echo "workflow set to " "${OPTARG}"
      [[ $_ngstype == 'GMS-ST' || $_ngstype == 'twist_ffpe_exome' || $_ngstype == 'test_panel' || $_ngstype == 'WGS' ]] || ( usage >&2; exit 1)
      ;;
    d)
      _analysis_dir=${OPTARG}
      echo "analysis dir set to " "${OPTARG}"
      ;;
    o)
      _out_dir=${OPTARG}
      echo "output dir set to " "${OPTARG}"
      ;;
     f)
      _tumor_fastq=${OPTARG}
      echo "tumor fastq is set to ${OPTARG}"
      ;;
     n)
      _normal_fastq=${OPTARG}
      echo "normal fastq is set ti ${OPTARG}"
      ;;
     r)
      rFlag=true;
      echo "Runflag is set to true"
      ;;
       *) echo "Invalid option: -${OPTARG}" >&2; usage >&2; exit 1;;
  esac
done

t_samplename=$(echo "$(basename ${_tumor_fastq})" | awk -F'_R_' '{print $1}')

# Function to create config
function balsamic_config() {
  balsamic config case \
    -t ${_tumor_fastq} \
    ${_normal_option} \
    --case-id ${_analysis}_${_ngstype}_${t_samplename} \
    --analysis-dir ${_analysis_dir} \
    -r ${_reference} \
    ${_panel_option} \
    --singularity ${_singularity} 
}

# Function to run balsamic analysis
function balsamic_run() {
  balsamic run analysis \
    -s ${_analysis_config} \
    -c ${_cluster_config} \
    --qos low \
    --profile qsub \
    --account batch.q ${_run_analysis} \
    --snakemake-opt '--latency-wait 30'
}

# Function to check completeness
function balsamic_endcheck() {
  balsamic report status \
    -s ${_analysis_config} \
    -m
}

module load miniconda/4.8.3
unset LD_PRELOAD
unset DISPLAY

if [[ ${_ngstype} == "twist_ffpe_exome" ]]; then
   _panel_option='-p /medstore/Development/GMS_solid_tumor/balsamic/beds/Twist_Exome_Target_hg19.bed'
elif [[ ${_ngstype} == "GMS-ST" ]]; then
   _panel_option='-p /absolute/path/to/panel/bed/here'
elif [[ ${_ngstype} == "test_panel" ]]; then 
   _panel_option='-p tests/test_data/references/panel/panel.bed'
else
  _panel_option=''
fi

if [[ ! -z ${_condaenv} ]]; then
  source activate ${_condaenv}
fi

if [[ -z ${_analysis_dir} ]]; then
  _analysis_dir='analysis_dir/'
  echo "analysis dir set to " "${_analysis_dir}"
fi

# Make sure _analysis_dir exists
mkdir -p ${_analysis_dir}

_prefix=/apps/bio/software/balsamic/BALSAMIC-6.0.0
_genome_ver=hg19
_cluster_config=${_prefix}/BALSAMIC/config/cluster.json
_singularity=${_prefix}/BALSAMIC/containers/balsamic_release_v6.0.0
#_reference=reference/${_genome_ver}/reference.json
_reference=${_prefix}/BALSAMIC_reference/hg19/reference.json
#_tumor_fastq=tests/test_data/fastq/S1_R_1.fastq.gz
#_normal_fastq=tests/test_data/fastq/S2_R_1.fastq.gz
_analysis_config=${_analysis_dir}'/'${_analysis}_${_ngstype}_${t_samplename}'/'${_analysis}_${_ngstype}_${t_samplename}'.json'

if [[ ! -z ${rFlag} ]]; then
  _run_analysis="-r"
fi

if [[ ${_analysis} == "TN" ]] && [[ -n "${_normal_fastq}" ]]; then
  _normal_option="-n ${_normal_fastq}"
else
  _normal_option=" "
fi

#    --snakemake-opt --singularity-args \
#    --snakemake-opt --cleanenv \

if [[ $_startmode == 'config' ]]; then
  balsamic_config
elif [[ $_startmode == 'run' ]]; then
  balsamic_run
else
  balsamic_config
  balsamic_run
fi
