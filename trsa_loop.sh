R1=''
R2=''
files=''

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
  printf "./trsa_loop -d <path/to/fastq> -core 4 -i <path/to/index> -o <output/file> -a <path/to/annotation>\n"
}

while getopts 'd:o:c:i:h' flag; do
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

for file in "$dir"/*.R1.fastq
do
  R1="$file"
  IFS="." read -r -a name <<< "$file"
  R2="${name[0]}"."${name[1]}".R2.fastq
  ./trsa.sh -s "$R1" -r "$R2" -c "$core" -i "$index"
done

# Generation de la table de correspondance entre ENCODE et HUGO
FILE=equivalence.txt;
if [ ! -d "$FILE" ]; then
  perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' \
    gencode.v24lift37.basic.annotation.gtf | sort | uniq > "$FILE"  ;
fi

# Construction de la table de comptage
# featureCounts

featureCounts -p -t exon -g gene_id -a "$annotation" \
              -o comptage.txt STAR/*.bam.bai

cat comptage.txt | sort > comptage_sort.txt

join -11 -11 equivalence.txt comptage_sort.txt > join.txt

awk '{print $2 "\t" $8"\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13}' join.txt > "$output"
