#!/usr/bin/env bash

# Note: can use 'snoderes -p PARTITION' to check resources but giving
#       a list is easier :)

#SBATCH --partition jic-long,jic-medium           # using SSD and only jic* nodes have them (spec. j512* nodes)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --time 1-0:00                             # time (D-H:m)
#SBATCH --output ./slurm_output/align.%A_%a.out   # STDOUT
#SBATCH --mem=40G                                 # memory pool for ALL cores
#SBATCH --array=1-10                              # default - override in 'sbatch' command
#SBATCH --gres=ssd
#SBATCH --localscratch=ssd:60                     # request 60GB SSD - enough?

# Minimise GC thread count (for fastqc and trimmomatic)
JAVA_OPTS="-XX:ParallelGCThreads=1"

cpus=${SLURM_CPUS_PER_TASK}
samthreads=6
ssd=${SLURM_LOCAL_SCRATCH}

function info() {
  echo "--> $(date "+%F %T") INFO : $@" >&2
}
function error() {
  echo "--> $(date "+%F %T") ERROR : $@" >&2
  exit 1
}

function resolve() {
  # return canonical file path

  local d
  local f
  d=$(dirname "$1")
  f=$(basename "$1")

  echo "$(cd "${d}" && pwd -P)/${f}"
}

function get_dtool_dataset_dir() {
  # return the directory of the dtool dataset (checks for ".dtool" directory)
  # returns "/" if file is not in a dataset

  local d=$1

  while [[ "${d}" != "/" ]]; do
    d=$(dirname ${d})

    if [[ -d ${d}/".dtool" ]]; then
      break
    fi
  done

  echo ${d}
}

echo "--------------------------------------------------------------------------------"

. "/hpc-home/${USER}/miniconda3/etc/profile.d/conda.sh"

CONDA_ENVS_PATH=/path/to/conda/envs

conda activate ${CONDA_ENVS_PATH}/conda-alignment-env

echo "software:"

echo "--> fastqc" && fastqc --version && echo
echo "--> hisat2" && hisat2 --version && echo
echo "--> samtools" && samtools --version && echo
echo "--> stringtie" && stringtie --version && echo
echo "--> trimmomatic" && trimmomatic -version && echo

echo "--------------------------------------------------------------------------------"

# enter bash "strict" mode
set -euo pipefail

# show trace
# set -x

fq_files=$1         # the file containing a list of fastq files/dirs
ref_fasta_file=$2
ref_gff_file=$3
ref_hisat2_index=$4
uniq_map_yn=$5

for f in "${fq_files}" "${ref_fasta_file}" "${ref_gff_file}"
do
  [[ ! -r "${f}" ]] && error "file is not readable (or it doesn't exist)! file = ${f}"
done

# get the fastq file/dir for this array job
sample_line=$(awk "NR == ${SLURM_ARRAY_TASK_ID}" ${fq_files})
sample_line=$(resolve $sample_line)

info "processing $fq_files (line ${SLURM_ARRAY_TASK_ID}) : $sample_line"

# are we in directory or file mode?
#
if [[ -d $sample_line ]]; then
  dir_mode=1
  sample_name=$(basename ${sample_line})

  # get the files in the sample directory
  samples=($(ls ${sample_line}/*_1.fq.gz))

  if [[ ${#samples[@]} -eq 0 ]]; then
    error "there are no files in the sample directory $sample_line"
  fi
elif [[ -f $sample_line ]]; then
  dir_mode=0
  sample_name=$(basename ${sample_line} | sed 's/_1.fq.gz//')

  # fastq files need to end with "_1.fq.gz"
  if [[ ${sample_line} != *"_1.fq.gz" ]]; then
    error "sample is invalid - it is expected to end with '_1.fq.gz'!"
  fi

  samples=("$sample_line")
else
  error "the sample line is not a file or directory : $sample_line"
fi

# require it to be in a 'dtool' dataset
dataset_dir=$(get_dtool_dataset_dir ${sample_line})
dataset_name=$(basename ${dataset_dir})

if [[ ${dataset_dir} == "/" ]]; then
  error "the sample is not part of a dtool dataset"
fi

echo "Arguments"
echo "========="
echo "sample name           : ${sample_name}"
echo "line from input file  : ${sample_line}"
echo "samples ['_1' only]   : ${samples[@]}"
echo "dataset name          : ${dataset_name} (dir=${dataset_dir})"
echo "reference fasta       : ${ref_fasta_file}"
echo "reference gff/gtf     : ${ref_gff_file}"
echo "hisat2 index          : ${ref_hisat2_index}"
echo "uniquely mapped reads : ${uniq_map_yn}"
echo " "

echo "--------------------------------------------------------------------------------"

info "processing..."

d=$(pwd -P)

qc_ds_dir="${d}/qc/${dataset_name}"

fastqc_raw_dir="${qc_ds_dir}/fastqc.raw"
fastqc_trm_dir="${qc_ds_dir}/fastqc.trim"
trimmo_log_dir="${qc_ds_dir}/trimmomatic"
hisat2_log_dir="${qc_ds_dir}/hisat2"
flagstat_dir="${qc_ds_dir}/flagstat"

intermediate_dir="${d}/intermediate_data/${dataset_name}"
final_dir="${d}/final_data/${dataset_name}"

mkdir -p "${fastqc_raw_dir}"
mkdir -p "${fastqc_trm_dir}"
mkdir -p "${trimmo_log_dir}"
mkdir -p "${hisat2_log_dir}"
mkdir -p "${flagstat_dir}"

mkdir -p "${intermediate_dir}"
mkdir -p "${final_dir}"

info "cleaning up old data for '${sample_name}' from intermediate directory..."

rm -f ${fastqc_raw_dir}/${sample_name:?}*   || true
rm -f ${fastqc_trm_dir}/${sample_name:?}*   || true
rm -f ${trimmo_log_dir}/${sample_name:?}*   || true
rm -f ${hisat2_log_dir}/${sample_name:?}*   || true
rm -f ${flagstat_dir}/${sample_name:?}*     || true

rm -f ${intermediate_dir}/${sample_name:?}* || true
rm -f ${final_dir}/${sample_name:?}*        || true

# go to the processing directory
#
cd $intermediate_dir

info "entering ssd directory ${ssd} ..."
cd $ssd

# for the trimmed fq files - will be passed to hisat2
samples_1=()
samples_2=()

# qc and trim each sample
#
for sample_1 in "${samples[@]}"
do
  echo "processing ${sample_1}"

  sample_id=$(basename ${sample_1} | sed 's/_1.fq.gz//')          # e.g. Rjap_1
  sample_2=$(echo ${sample_1} | sed 's/_1.fq.gz/_2.fq.gz/')       # e.g. /path/to/Rjap_1_2.fq.gz

  echo "- sample name     : ${sample_name}" 
  echo "- sample id       : ${sample_id}" 
  echo "- sample 1        : ${sample_1}"
  echo "- sample 2        : ${sample_2}"

  info "running fastqc - initial raw read QC..."

  srun fastqc \
    -t ${cpus} \
    -o ${fastqc_raw_dir} \
    ${sample_1} \
    ${sample_2}

  info "trimming with Trimmomatic..."

  trimmo_adapters=${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE-2.fa

  trimmo_opts="ILLUMINACLIP:${trimmo_adapters}:2:30:10:1:true "
  trimmo_opts+="HEADCROP:15 "
  trimmo_opts+="SLIDINGWINDOW:4:15 "
  trimmo_opts+="MINLEN:50 "

  # P files - paired reads after trimming
  # U files - unpaired reads after trimming

  trim_1p=${sample_id}.trimmed_1P.fq
  trim_1u=${sample_id}.trimmed_1U.fq
  trim_2p=${sample_id}.trimmed_2P.fq
  trim_2u=${sample_id}.trimmed_2U.fq

  srun trimmomatic \
    PE \
    -threads ${cpus} \
    -phred33 \
    ${sample_1} \
    ${sample_2} \
    ${trim_1p} \
    ${trim_1u} \
    ${trim_2p} \
    ${trim_2u} \
    ${trimmo_opts} 2> ${trimmo_log_dir}/${sample_id}.log

  info "zipping trimmed fastq files..."

  [[ -f "${trim_1p}" ]] && srun pigz -p ${cpus} -f ${trim_1p}
  [[ -f "${trim_2p}" ]] && srun pigz -p ${cpus} -f ${trim_2p}

  trim_1p=${trim_1p}.gz
  trim_2p=${trim_2p}.gz

  samples_1+=("${trim_1p}")
  samples_2+=("${trim_2p}")

  info "running fastqc - trimmed reads..."

  srun fastqc \
    -t ${cpus} \
    -o ${fastqc_trm_dir} \
    ${trim_1p} \
    ${trim_2p}
done

# Parameters for HISAT2 from:
#   Baruzzo et al. Simulation-based comprehensive benchmarking of RNA-seq aligners. Nat Methods 14, 135–139 (2017).
#   https://doi.org/10.1038/nmeth.4106
#   Supplementary Note 3
#
# --pen-noncansplice <int>    penalty for a non-canonical splice site (12)
# --mp <int>,<int>            max and min penalties for mismatch; lower qual = lower penalty <6,2>

info "aligning with hisat2..."

extra_options="--pen-noncansplice 20 --mp 1,0 "

# align multiple files using "-1 file1_1,file2_1 -2 file1_2,file2_2" so..
# ..convert array to comma-delimited variables
trimmed_1=$(IFS=, ; echo "${samples_1[*]}")
trimmed_2=$(IFS=, ; echo "${samples_2[*]}")

bam_aligned=${sample_name}.aligned.bam

srun hisat2 \
  --threads ${cpus} \
  --new-summary \
  --summary-file ${hisat2_log_dir}/${sample_name}.log \
  --time \
  --dta \
  --rna-strandness RF \
  ${extra_options} \
  -x ${ref_hisat2_index} \
  -1 ${trimmed_1} \
  -2 ${trimmed_2} | samtools sort -@ ${samthreads} -o ${bam_aligned} -O bam -
  
conda deactivate

# For full details of SAM and BAM format see:
#   https://samtools.github.io/hts-specs/SAMv1.pdf
#
# To decode SAM flags see:
#   https://broadinstitute.github.io/picard/explain-flags.html
#   https://www.samformat.info/sam-format-flag
#
# For example:
#
#  flag   meaning
#  ====   =======
#     1   Read paired
#     2   Read mapped in proper pair ('proper' refers to read orientation)
#     4   Read unmapped
#     8   Mate unmapped
#    12   Read + mate unmapped (=4+8)
#   256   Not primary alignment
#   260   Read unmapped + not primary alignment (=256 + 4)
#   264   Mate unmapped + not primary alignment (=256 + 8)
#
# Note: primary alignment => only one primary allowed
#       supplemental alignment => something to do with chimeric alignments
#
# When using samtools
#   -f INT   only include reads with all  of the FLAGs in INT present [0]
#   -F INT   only include reads with none of the FLAGS in INT present [0]

# use more recent samtools - more subcommands have 'samthreads' plus the ability to filter on tags
conda activate ${CONDA_ENVS_PATH}/conda-samtools-env

info "indexing..."
srun samtools index -@ ${samthreads} ${bam_aligned}

info "generating flagstat alignment stats for QC (sorted BAM)..."
flagstat_file=${flagstat_dir}/${sample_name}.align.log
srun samtools flagstat -@ ${samthreads} ${bam_aligned} > ${flagstat_file}

info "filtering aligned reads (keep proper pairs)..."
bam_mapped=${sample_name}.mapped-multi.bam
srun samtools view -@ ${samthreads} -bh -f 2 -O bam ${bam_aligned} | samtools sort -@ ${samthreads} -o ${bam_mapped} -

info "indexing..."
srun samtools index -@ ${samthreads} ${bam_mapped}

info "running flagstat for mapped reads (some multimapped)..."
flagstat_file=${flagstat_dir}/${sample_name}.multi.log
srun samtools flagstat -@ ${samthreads} ${bam_mapped} > ${flagstat_file}

info "extracting uniquely mapped reads..."
bam_mapped_unique=${sample_name}.mapped-unique.bam
srun samtools view -@ ${samthreads} -h --tag NH:1 ${bam_mapped} -o ${bam_mapped_unique}

info "running flagstat for uniquely mapped reads..."
flagstat_file=${flagstat_dir}/${sample_name}.unique.log
srun samtools flagstat -@ ${samthreads} ${bam_mapped_unique} > ${flagstat_file}

# Extract unmapped reads
#
# More details: http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/1040-2/
#
info "extracting unmapped reads..."

bam_unmapped=${sample_name}.unmapped.bam

# an unmapped read whose mate is mapped
srun samtools view -@ ${samthreads} -u -f 4 -F 264 ${bam_aligned} -o ${bam_unmapped}.1.bam

# a mapped read who's mate is unmapped
srun samtools view -@ ${samthreads} -u -f 8 -F 260 ${bam_aligned} -o ${bam_unmapped}.2.bam

# both reads of the pair are unmapped
srun samtools view -@ ${samthreads} -u -f 12 -F 256 ${bam_aligned} -o ${bam_unmapped}.3.bam

srun samtools merge -@ ${samthreads} -u - ${bam_unmapped}.[123].bam | samtools sort -@ ${samthreads} -o ${bam_unmapped} -

conda deactivate

info "copying bam files to the intermediate directory..."

srun cp $bam_mapped         $intermediate_dir/
srun cp $bam_mapped_unique  $intermediate_dir/
srun cp $bam_unmapped       $intermediate_dir/

conda activate ${CONDA_ENVS_PATH}/conda-alignment-env

# Assemble & quantify reads
#
# -A gene abundance estimation output file
# -B enable output of Ballgown table files which will be created in the same directory as the
#    output GTF (requires -G, -o recommended)
# -e only estimate the abundance of given reference transcripts (requires -G)
# -G reference annotation to use for guiding the assembly process (GTF/GFF3)
# -l name prefix for output transcripts
# -o output path/file name for the assembled transcripts GTF
# -p number of threads (CPUs) to use
#
# Note: stringtie attempts to correct for multi-mapped reads by assigning proportionately to the transcript (turn off with "-u")
#
# -u  Turn off multi-mapping correction. In the default case this correction is enabled, and each read
#     that is mapped in n places only contributes 1/n to the transcript coverage instead of 1. 
#
info "assemble and quantify with stringtie..."

ds_sample=${dataset_name}_${sample_name}

stringtie_outdir=${final_dir}/ballgown/ERR_${ds_sample}
stringtie_genes=${stringtie_outdir}/${ds_sample}_gene_abundance.tab
stringtie_transcripts=${stringtie_outdir}/${ds_sample}_transcriptome.gtf

if [[ ${uniq_map_yn} == "y" ]]; then
  bam=${bam_mapped_unique}
else
  bam=${bam_mapped}
fi

srun stringtie \
  ${bam} \
  -A ${stringtie_genes} \
  -B \
  -e \
  -G ${ref_gff_file} \
  -l ${sample_name} \
  -o ${stringtie_transcripts} \
  -p ${cpus}

info "successfully produced ${stringtie_transcripts}"

echo "--------------------------------------------------------------------------------"

info "success - aligned ${sample_name}"

echo "--------------------------------------------------------------------------------"
