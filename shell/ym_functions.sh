##########################################################################@ # 合并多个RNA-seq的count文件
function func_CombineRNAseqCountFilesFromSingleSample(){
    if [ $# -eq 0 ];then
        echo "check arguments!"
        exit 1
    fi

    while getopts "d:c:m:l:" arg; do
        case "$arg" in
            d)
                inputCountDir="$OPTARG"
                ;;
            c)
                countsFilename="$OPTARG"
                ;;
            m)
                metainfoFilename="$OPTARG"
                ;;
            l)
                librarys="$OPTARG"
                ;;
            ?)
                echo "Invalid option: -$OPTARG"
                exit 1
                ;;
        esac
    done

    echo -e "inputCountDir: $inputCountDir\n\nlibrarys: $librarys\n"
    echo -e "output: \n\t$countsFilename\n\t$metainfoFilename"
    
    for i in $librarys;do
        awk 'NR==1{print "ENSG\t"filename; print $0}NR>1{print $0}' filename=$i ${inputCountDir}/${i}/${i}.count >tmp.$i.count
    done
    # 生成counts文件
    paste -d "\t" tmp.*.count | awk '{
        for(i=1;i<=NF;i++){
            if(i==1) printf $i"\t"
            if(i%2==0 && i<NF) printf $i"\t"
            if(i%2==0 && i==NF) printf $i}
        printf "\n"
        }' > ${countsFilename}
    # 生成metainfo文件
    ls tmp.*.count | awk '{
        gsub(/tmp.|.count/, "")
        if(NR==1) print "sampleLibrary\tsample\tcondition\n"$0
        if(NR>1) print $0
        }' > ${metainfoFilename}
    rm tmp.*.count

    echo -e "\n......根据情况补全metainfo文件\n格式（制表符分隔）：文库名称 文库短名称 文库组别名称"
}













##########################################################################@ # rename Dir and files inner #####
function func_rename_Dir_Files_inside(){
    oldName=$1
    newName=$2

    mv $oldName $newName
    cd $newName
    find . -name "*$oldName*" | sed -e "p;s/$oldName/$newName/" | xargs -n2 mv
    cd ..
}
### usage: (in target dir) func_rename_Dir_Files_inside $old $new












##########################################################################@ # rpkm files #####
function func_normalize_by_rpkm(){
    inputNodupBam=$1
    species=$2

    if [[ "${species}" == "hg38" ]]; then
        blackListFile="/data2/yuanming/refGenome/annotation/GRCh38_unified_blacklist.bed"
    elif [[ "${species}" == "mm10" ]]; then
        blackListFile="/data2/yuanming/refGenome/annotation/mm10-blacklist.v2.bed"
    fi

    if [[ ! -f "${inputNodupBam}.bai" ]]; then
        echo "Index : ${inputNodupBam}"
        samtools index -@ 20 ${inputNodupBam};
    fi

    echo -e "\n\n\nNormalizing : ${inputNodupBam/.nodup.bam/}"    
    /home/yuanming/mambaforge/bin/bamCoverage \
    --extendReads 150 -p 20 \
    --blackListFileName ${blackListFile} \
    --binSize 10 \
    --normalizeUsing RPKM \
    -b $inputNodupBam \
    -o ${inputNodupBam/.bam/}.RPKM.bigwig
}
### usage: func_normalize_by_rpkm $inputNodupBam $species











##########################################################################@ # featureCount #####
function fun_featureCount_multiBigWig(){
    gtfFile=$1
    assayName=$2
    PREFIXsample=$3
    inputDir=$4

    echo "make sure "
    featureCounts -T 20 -p -t peak -g gene_id -a $gtfFile -o $assayName.count $(for i in $PREFIXsample;do echo $inputDir/$i/$i.nodup.bam;done | xargs)
    echo -e "Geneid\t$(echo $PREFIXsample | sed -E 's/\s/\t/g')" > ${assayName}.de_counts.txt
    awk 'NR>2{ for (i=1; i<=NF; i++) if (i < 2 || i > 6) printf "%s%s", $i, (i == NF ? "\n" : "\t") }' $assayName.count >> ${assayName}.de_counts.txt
    echo -e "sampleLibrary\tsample\tcondition\n$(echo $PREFIXsample | sed -E 's/\s/\n/g')" >${assayName}.de_metainfo.txt
    echo -e "\n...complete metainfo manually...\n"
}
### usage: fun_featureCount_multiBigWig "$gtfFile" "$assayName" "$PREFIXsample" "$inputDir"










# **************************************
# BED file operation
# **************************************
##########################################################################@ # transform to bigbed #####
function fun_trans_to_bigBed(){
    inputBed=$1
    species=$2

    if [[ "${species}" == "hg38" ]]; then
        GENOME_SIZE=/data/public/refGenome/bwa_index/hg38/ChromInfo.txt
    elif [[ "${species}" == "mm10" ]]; then
        GENOME_SIZE=/data/public/refGenome/bwa_index/mm10/mm10.chrom.sizes
    fi

    sort -k1,1 -k2,2n ${inputBed} > ${inputBed/.bed/.sorted.bed}
    bedToBigBed ${inputBed/.bed/.sorted.bed} ${GENOME_SIZE} ${inputBed/.bed/.bigbed}

    rm ${inputBed/.bed/.sorted.bed}
    echo "Done: ${inputBed} ---> ${inputBed/.bed/.bigbed}"
}
### usage: fun_trans_to_bigBed inputBed.bed hg38


function fun_index_bed(){
    inputBed=$1

    sort -k1,1 -k2,2n ${inputBed} > ${inputBed/.bed/.sorted.bed}
    bgzip ${inputBed/.bed/.sorted.bed}
    tabix -p bed ${inputBed/.bed/.sorted.bed}.gz
}
### usage: fun_index_bed inputBed.bed









# *******************************************
# transform bedgraph to bigwig
# *******************************************
fun_trans_bedGraph_to_Bigwig(){
    input_bg=$1
    REF=$2

    if [ "${REF}" == "hg38" ];
    then
        GENOME_SIZE=/data/public/refGenome/bwa_index/hg38/ChromInfo.txt

    elif [ "${REF}" == "mm10" ];
    then
        GENOME_SIZE=/data/public/refGenome/bwa_index/mm10/mm10.chrom.sizes

    else
        echo -e "Received an unexpected reference genome: ${REF}"
        exit
    fi

    echo -e "\n\n\nTrans $input_bg"
    bedGraphToBigWig ${input_bg} ${GENOME_SIZE} ${input_bg/.bedgraph/.bigwig}
}
### usage: fun_trans_bedGraph_to_Bigwig $input_bg $REF
