R1=''
R2=''
files=''

print_usage() {
  printf "Usage: ...\n"
}

while getopts 's:r:f:c:i:a:h' flag; do
  case "${flag}" in
    s) R1="${OPTARG}" ;;
    r) R2="${OPTARG}" ;;
    c) core="${OPTARG}" ;;
    i) index="${OPTARG}" ;;
    f) files="${OPTARG}" ;;
    a) annotation="${OPTARG}" ;;
    h) print_usage
       exit 1 ;;
  esac
done

echo $R1
echo $R2;
echo $files;

DIR=trimmomatic;
if [ ! -d "$DIR" ]; then
  mkdir $DIR;
fi

DIR=STAR;
if [ ! -d "$DIR" ]; then
  mkdir $DIR;
fi

DIR=featureCounts;
if [ ! -d "$DIR" ]; then
  mkdir $DIR;
fi

# Activation de conda
conda activate

# Analyse fastqc
fastqc "$R1" "$R2"

# Trimming
IFS="/" read -r -a name <<< "$R1"
R1_name=${name[-1]}
IFS='.' read -r -a array <<< "$R1_name"
basename=${array[0]}
trimmomatic PE "$R1" "$R2" -baseout trimmomatic/"$basename".fastq LEADING:20 TRAILING:20 MINLEN:5
echo -----------------------------------------
echo $basename
echo -----------------------------------------
# STAR
STAR --runThreadN "$core" --outFilterMultimapNmax 1 --genomeDir "$index" \
     --outSAMattributes All --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix STAR/"$basename" \
     --readFilesIn trimmomatic/"$basename"_1P.fastq trimmomatic/"$basename"_2P.fastq ;

# samtools
echo samtools pour "$basename"
samtools index STAR/"$basename"*.bam

