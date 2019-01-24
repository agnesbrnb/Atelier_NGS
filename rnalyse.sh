
print_usage() {
  printf "Programme permettant de faire une analyse de données RNAseq\n"
  printf "Usage :\n"
  printf "  [-d] chemin vers le repertoire contenant les fichiers fastq\n"
  printf "  [-c] nombre de coeur(s) utilisés\n"
  printf "  [-i] chemin vers le répertoire contenant les index du génome\n"
  printf "  [-o] nom du fichier output\n"
  printf "  [-a] chemin vers les annotations chromosome"
  printf "\n"
  printf "Exemple:\n"
  printf "./rnanalyse -d <path/to/fastq> -core 4 -i <path/to/index> -o <output/file> -a <path/to/annotation>\n"
}

while getopts 'd:o:c:i:a:h' flag; do
  case "${flag}" in
    d) dir="${OPTARG}" ;;
    c) core="${OPTARG}" ;;
    i) index="${OPTARG}" ;;
    o) output="${OPTARG}" ;;
    a) annotation="${OPTARG}" ;;
    h) print_usage
       exit 1 ;;
  esac
done

# Création des repertoires pour les fastq "trimmer"

DIR=trimmomatic;
if [ ! -d "$DIR" ]; then
  mkdir $DIR;
fi

# Création deu repertoires contenant les fichiers BAM et BAM indexer - BAI

DIR=STAR;
if [ ! -d "$DIR" ]; then
  mkdir $DIR;
fi

for file in "$dir"/*.R1.fastq
do
  R1="$file"
  IFS="." read -r -a name <<< "$file"
  R2="${R1//R1/R2}"


  IFS="/" read -r -a name <<< "$R1"
  R1_name=${name[-1]}
  IFS='.' read -r -a array <<< "$R1_name"
  basename=${array[0]}

  # Analyse fastqc
  fastqc "$R1" "$R2"

  # Trimming
  trimmomatic PE "$R1" "$R2" -baseout trimmomatic/"$basename".fastq LEADING:20 TRAILING:20 MINLEN:5
  
  # Mapping des sequences avec la commande STAR
  # contient les fichiers BAM
  STAR --runThreadN "$core" --outFilterMultimapNmax 1 --genomeDir "$index" \
       --outSAMattributes All --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix STAR/"$basename" \
       --readFilesIn trimmomatic/"$basename"_1P.fastq trimmomatic/"$basename"_2P.fastq ;

  # Indexation des fichiers BAM avec samtools
  samtools index STAR/"$basename"*.bam

done

# Generation de la table de correspondance entre ENCODE et HUGO
FILE=equivalence.txt;
if [ ! -f "$FILE" ]; then
  perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' \
    "$annotation" | sort | uniq > "$FILE"  ;
fi

# Construction de la table de comptage
# featureCounts

featureCounts -p -t exon -g gene_id -a "$annotation" \
              -o comptage.txt STAR/*.bam

cat comptage.txt | sort > comptage_sort.txt

join equivalence.txt comptage_sort.txt > join.txt

awk '{print $2 "\t" $8"\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13}' join.txt > "$output"

