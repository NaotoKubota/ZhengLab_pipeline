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
 ____        __           __                                           
/\  _`\     /\ \__       /\ \                                          
\ \ \L\_\ __\ \ ,_\   ___\ \ \___                                      
 \ \  _\/'__`\ \ \/  /'___\ \  _ `\                                    
  \ \ \/\  __/\ \ \_/\ \__/\ \ \ \ \                                   
   \ \_\ \____\\ \__\ \____\\ \_\ \_\                                  
    \/_/\/____/ \/__/\/____/ \/_/\/_/                                  
EOF

echo "
Version: ${VERSION}
"

function usage {
    cat << EOS

ZhengLab_Fetch ${VERSION} -Pipeline for fetching NGS data in ZhengLab-

Usage: $(basename "$0") -i ID -o output_dir

    -h  Display help
    -v  Show version
    -i  ID
    -o  Output directory

EOS
    exit 2
}

# version
function version {
    cat << EOS

ZhengLab_Fetch ${VERSION} -Pipeline for fetching NGS data in ZhengLab-

EOS
    exit 2
}

function lack_of_necessary_param() {
    usage
    exit 1
}

IS_THERE_NECESSARY_OPT_i=false
IS_THERE_NECESSARY_OPT_o=false

while getopts "i:o:vh" optKey; do
    case "$optKey" in
    i)
        IS_THERE_NECESSARY_OPT_i=true
        ID=${OPTARG}
        ;;
    o)
        IS_THERE_NECESSARY_OPT_o=true
        OUTPUTDIR=${OPTARG}
        ;;
    v|version)
        version
        ;;
    h|* )
        usage
        ;;
    esac
done


if [ "${IS_THERE_NECESSARY_OPT_i}" == false ] || [ "${IS_THERE_NECESSARY_OPT_o}" == false ]; then
    lack_of_necessary_param
fi;


# print parameters
echo -e "
ID:\t${ID}
OUTPUTDIR:\t${OUTPUTDIR}
"


# make output directory
mkdir -p ${OUTPUTDIR}/fastq/log


# fetch metadata
echo -e "Fetching metadata ..."
# check if file exists
if [ ! -f ${OUTPUTDIR}/${ID}_metadata.json ]; then

    echo -e "Downloading ${ID}_metadata.json ..."
    ffq ${ID} -o ${OUTPUTDIR}/${ID}_metadata.json

else

    echo -e "${ID}_metadata.json already exists!"

fi


# fetch fastq files
echo -e "Fetching fastq files ..."
# check if file exists
if [ ! -f ${OUTPUTDIR}/${ID}_fastq_ftp.txt ]; then

    echo -e "Downloading ${ID}_fastq_ftp.txt ..."
    ffq --ftp ${ID} | jq '.[] | .md5 + " " + .url' -r | grep -e "fastq.gz$" > ${OUTPUTDIR}/${ID}_fastq_ftp.txt

else

    echo -e "${ID}_fastq_ftp.txt already exists!"

fi

rm -rf ${OUTPUTDIR}/fastq/log/md5sum.log && \
touch ${OUTPUTDIR}/fastq/log/md5sum.log
cat ${OUTPUTDIR}/${ID}_fastq_ftp.txt | while read line
do

	# filename
	filename=$(echo -e ${line} | cut -d" " -f 2 | xargs -n1 basename)
    # URL
    URL=$(echo -e ${line} | cut -d" " -f 2)
    # md5
    md5=$(echo -e ${line} | cut -d" " -f 1)

	# check if file exists
	if [ ! -f ${OUTPUTDIR}/fastq/${filename} ]; then

		echo -e "Downloading ${filename} ..."
		wget -P ${OUTPUTDIR}/fastq/ ${URL} 2>> ${OUTPUTDIR}/fastq/log/${filename}.log && \
        # md5
        echo -e "${md5} ${OUTPUTDIR}/fastq/${filename}" | md5sum -c - >> ${OUTPUTDIR}/fastq/log/md5sum.log

	else

		echo -e "${filename} already exists!"

	fi

done


echo -e \
"
Finished!
OUTPUTDIR:\t${OUTPUTDIR}
"
