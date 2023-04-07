#!/bin/bash -e
set -o pipefail
# Prevent commands misbehaving due to locale differences
export LC_ALL=C


VERSION="v0.1"

cat << "EOF"

  ________   __                             __                __        
/\_____  \ /\ \                           /\ \              /\ \       
\/____//'/'\ \ \___      __    ___      __\ \ \         __  \ \ \____  
     //'/'  \ \  _ `\  /'__`\/' _ `\  /'_ `\ \ \  __  /'__`\ \ \ '__`\ 
    //'/'___ \ \ \ \ \/\  __//\ \/\ \/\ \L\ \ \ \L\ \/\ \L\.\_\ \ \L\ \
    /\_______\\ \_\ \_\ \____\ \_\ \_\ \____ \ \____/\ \__/.\_\\ \_,__/
    \/_______/ \/_/\/_/\/____/\/_/\/_/\/___L\ \/___/  \/__/\/_/ \/___/ 
                                        /\____/                        
                                        \_/__/                         
 ____    __  __  ______                                                
/\  _`\ /\ \/\ \/\  _  \                                               
\ \ \L\ \ \ `\\ \ \ \L\ \             ____     __     __               
 \ \ ,  /\ \ , ` \ \  __ \  _______  /',__\  /'__`\ /'__`\             
  \ \ \\ \\ \ \`\ \ \ \/\ \/\______\/\__, `\/\  __//\ \L\ \            
   \ \_\ \_\ \_\ \_\ \_\ \_\/______/\/\____/\ \____\ \___, \           
    \/_/\/ /\/_/\/_/\/_/\/_/         \/___/  \/____/\/___/\ \          
                                                         \ \_\         
                                                          \/_/         
EOF

echo "Version: ${VERSION}"

function usage {
    cat << EOS

ZhengLab_RNAseq ${VERSION} -Pipeline for RNA-seq analysis in ZhengLab-

Usage: $(basename "$0") -i fastq.txt -c config.txt

    -h  Display help
    -v  Show version
    -i  Fastq path table
    -c  Config file

EOS
    exit 2
}

# version
function version {
    cat << EOS

ZhengLab_RNAseq ${VERSION} -Pipeline for RNA-seq analysis in ZhengLab-

EOS
    exit 2
}

function lack_of_necessary_param() {
    usage
    exit 1
}

IS_THERE_NECESSARY_OPT_i=false
IS_THERE_NECESSARY_OPT_c=false

while getopts "i:c:vh" optKey; do
    case "$optKey" in
    i)
        IS_THERE_NECESSARY_OPT_i=true
        FASTQ_TABLE=${OPTARG}
        ;;
    c)
        IS_THERE_NECESSARY_OPT_c=true
        CONFIG=${OPTARG}
        ;;
    v|version)
        version
        ;;
    h|* )
        usage
        ;;
    esac
done


if [ "${IS_THERE_NECESSARY_OPT_i}" == false ] || [ "${IS_THERE_NECESSARY_OPT_c}" == false ]; then
    lack_of_necessary_param
fi;


# source config file
source ${CONFIG}


# Check if the necessary parameters are set
if [ -z "${NUM_PROCESS}" ]; then
	echo "NUM_PROCESS is not set. Set NUM_PROCESS in config file."
	exit 1
fi

if [ -z "${STAR_INDEX}" ]; then
	echo "STAR_INDEX is not set. Set STAR_INDEX in config file."
	exit 1
fi

if [ -z "${DEEPTOOLS_STRANDED}" ]; then
	echo "DEEPTOOLS_STRANDED is not set. Set DEEPTOOLS_STRANDED in config file."
	exit 1
fi

if [ -z "${OUTPUT}" ]; then
	echo "OUTPUT is not set. Set OUTPUT in config file."
	exit 1
fi

if [ -z "${SKIP_STEP1}" ]; then
	echo "SKIP_STEP1 is not set. Set SKIP_STEP1 in config file."
	exit 1
fi

if [ -z "${SKIP_STEP2}" ]; then
	echo "SKIP_STEP2 is not set. Set SKIP_STEP2 in config file."
	exit 1
fi

if [ -z "${SKIP_STEP3}" ]; then
	echo "SKIP_STEP3 is not set. Set SKIP_STEP3 in config file."
	exit 1
fi


# print parameters
echo -e "
Input:\t${FASTQ_TABLE}
Config:\t${CONFIG}
CPU:\t${NUM_PROCESS}
Output:\t${OUTPUT}
"


# execution time
SECONDS=0


# make output directory
mkdir -p ${OUTPUT}/MultiQC


# make log directory
mkdir -p ${OUTPUT}/log/MultiQC


## Step1: fastp
if [ "${SKIP_STEP1}" == false ]; then

	echo -e "[fastp] Started."

	# make output directory
	echo -e "[fastp] Make output directory: ${OUTPUT}/fastp"
	mkdir -p ${OUTPUT}/fastp

	# make log directory
	echo -e "[fastp] Make log directory: ${OUTPUT}/log/fastp"
	mkdir -p ${OUTPUT}/log/fastp


	cat ${FASTQ_TABLE} | while read line || [ -n "${line}" ]
	do

		READTYPE=$(echo -e "${line}" | cut -f 1)
		ID=$(echo -e "${line}" | cut -f 2)
		echo -e "[fastp] Run ${ID}.."

		if [ "${READTYPE}" == "paired" ]; then

			R1=$(echo -e "${line}" | cut -f 3)
			R2=$(echo -e "${line}" | cut -f 4)

			fastp \
			-i ${R1} \
			-I ${R2} \
			-o ${OUTPUT}/fastp/${ID}_1_filtered.fastq.gz \
			-O ${OUTPUT}/fastp/${ID}_2_filtered.fastq.gz \
			-h ${OUTPUT}/log/fastp/${ID}_fastp.html \
			-j ${OUTPUT}/log/fastp/${ID}_fastp.json \
			${FASTP_ARGS} \
			-w ${NUM_PROCESS} 2>> ${OUTPUT}/log/fastp/fastp.log

		elif [ "${READTYPE}" == "single" ]; then

			R1=$(echo -e "${line}" | cut -f 3)

			fastp \
			-i ${R1} \
			-o ${OUTPUT}/fastp/${ID}_filtered.fastq.gz \
			-h ${OUTPUT}/log/fastp/${ID}_fastp.html \
			-j ${OUTPUT}/log/fastp/${ID}_fastp.json \
			${FASTP_ARGS} \
			-w ${NUM_PROCESS} 2>> ${OUTPUT}/log/fastp/fastp.log

		else

			echo -e "Invalid read type: ${READTYPE}"
			exit 1

		fi

	done

	# Copy log files for MultiQC
	cp ${OUTPUT}/log/fastp/*_fastp.json ${OUTPUT}/MultiQC/

else

    echo -e "[fastp] Skipped."

fi

i_one=${SECONDS}
sec=$((i_one % 60))
min=$(((i_one % 3600) / 60))
hrs=$((i_one / 3600))
timestamp=$(printf "%d:%02d:%02d" "$hrs" "$min" "$sec")
echo -e "[fastp] Execution time: ${timestamp}"

## Step2: STAR
if [ "${SKIP_STEP2}" == false ]; then

	echo -e "[STAR] Started."

	# make output directory
	echo -e "[STAR] Make output directory: ${OUTPUT}/STAR"
	mkdir -p ${OUTPUT}/STAR

	# make log directory
	echo -e "[STAR] Make log directory: ${OUTPUT}/log/STAR"
	mkdir -p ${OUTPUT}/log/STAR

	cat ${FASTQ_TABLE} | while read line || [ -n "${line}" ]
	do

		READTYPE=$(echo -e "${line}" | cut -f 1)
		ID=$(echo -e "${line}" | cut -f 2)
		echo -e "[STAR] Run ${ID}.."
		mkdir -p ${OUTPUT}/STAR/${ID}

		if [ "${READTYPE}" == "paired" ]; then

			STAR \
            --runThreadN ${NUM_PROCESS} \
            --genomeDir ${STAR_INDEX} \
            --readFilesIn \
            ${OUTPUT}/fastp/${ID}_1_filtered.fastq.gz \
            ${OUTPUT}/fastp/${ID}_2_filtered.fastq.gz \
            --readFilesCommand zcat \
            ${STAR_ARGS} \
            --outFileNamePrefix \
            ${OUTPUT}/STAR/${ID}/${ID}_ 2>> ${OUTPUT}/log/STAR/STAR.log

		elif [ "${READTYPE}" == "single" ]; then

			STAR \
            --runThreadN ${NUM_PROCESS} \
            --genomeDir ${STAR_INDEX} \
            --readFilesIn \
            ${OUTPUT}/fastp/${ID}_filtered.fastq.gz \
            --readFilesCommand zcat \
            ${STAR_ARGS} \
            --outFileNamePrefix \
            ${OUTPUT}/STAR/${ID}/${ID}_ 2>> ${OUTPUT}/log/STAR/STAR.log


		else

			echo -e "Invalid read type: ${READTYPE}"
			exit 1

		fi

		echo -e "[Samtools] Convert SAM to BAM: ${ID}"
		samtools sort \
		-@ ${NUM_PROCESS} \
		-O bam \
		-o ${OUTPUT}/STAR/${ID}/${ID}.sort.bam \
		${OUTPUT}/STAR/${ID}/${ID}_Aligned.out.sam && \
		samtools index ${OUTPUT}/STAR/${ID}/${ID}.sort.bam && \
		rm -rf ${OUTPUT}/STAR/${ID}/${ID}_Aligned.out.sam

	done

	# Copy log files for MultiQC
	cp ${OUTPUT}/STAR/*/*_Log.final.out ${OUTPUT}/MultiQC/

else

	echo -e "[STAR] Skipped."

fi

i_two=$((${SECONDS} - ${i_one}))
sec=$((i_two % 60))
min=$(((i_two % 3600) / 60))
hrs=$((i_two / 3600))
timestamp=$(printf "%d:%02d:%02d" "$hrs" "$min" "$sec")
echo -e "[STAR] Execution time: ${timestamp}"


## Step3: deepTools bamCoverage
if [ "${SKIP_STEP3}" == false ]; then

	echo -e "[deepTools bamCoverage] Started."

	# make output directory
	echo -e "[deepTools bamCoverage] Make output directory: ${OUTPUT}/deepTools"
	mkdir -p ${OUTPUT}/deepTools

	# make log directory
	echo -e "[deepTools bamCoverage] Make log directory: ${OUTPUT}/log/deepTools"
	mkdir -p ${OUTPUT}/log/deepTools

	cat ${FASTQ_TABLE} | while read line || [ -n "${line}" ]
	do

		READTYPE=$(echo -e "${line}" | cut -f 1)
		ID=$(echo -e "${line}" | cut -f 2)
		echo -e "[deepTools bamCoverage] Run ${ID}.."

		if [ ${DEEPTOOLS_STRANDED} == true ]; then

			echo -e "[deepTools bamCoverage] Forward strand.."
			bamCoverage \
			-p ${NUM_PROCESS} \
			${DEEPTOOLS_ARGS} \
			--filterRNAstrand forward \
			-b ${OUTPUT}/STAR/${ID}/${ID}.sort.bam \
			-o ${OUTPUT}/deepTools/${ID}_forward.bw \
			2>> ${OUTPUT}/log/deepTools/deepTools.log

			echo -e "[deepTools bamCoverage] Reverse strand.."
			bamCoverage \
			-p ${NUM_PROCESS} \
			${DEEPTOOLS_ARGS} \
			--filterRNAstrand reverse \
			-b ${OUTPUT}/STAR/${ID}/${ID}.sort.bam \
			-o ${OUTPUT}/deepTools/${ID}_reverse.bw \
			2>> ${OUTPUT}/log/deepTools/deepTools.log

		else

			echo -e "[deepTools bamCoverage] Unstranded.."
			bamCoverage \
			-p ${NUM_PROCESS} \
			${DEEPTOOLS_ARGS} \
			-b ${OUTPUT}/STAR/${ID}/${ID}.sort.bam \
			-o ${OUTPUT}/deepTools/${ID}.bw \
			2>> ${OUTPUT}/log/deepTools/deepTools.log

		fi

	done

else

	echo -e "[deepTools bamCoverage] Skipped."

fi

i_three=$((${SECONDS} - ${i_two} - ${i_one}))
sec=$((i_three % 60))
min=$(((i_three % 3600) / 60))
hrs=$((i_three / 3600))
timestamp=$(printf "%d:%02d:%02d" "$hrs" "$min" "$sec")
echo -e "[deepTools bamCoverage] Execution time: ${timestamp}"


# MultiQC
echo -e "[MultiQC] Started."
multiqc \
-o ${OUTPUT}/MultiQC \
${OUTPUT}/MultiQC/ \
2>> ${OUTPUT}/log/MultiQC/MultiQC.log


i=${SECONDS}
sec=$((i % 60))
min=$(((i % 3600) / 60))
hrs=$((i / 3600))
timestamp=$(printf "%d:%02d:%02d" "$hrs" "$min" "$sec")

echo -e \
"
Finished!
Total execution time: ${timestamp}
"
