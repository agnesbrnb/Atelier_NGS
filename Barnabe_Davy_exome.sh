#bwa index -a bwtsw chr16.fa.gz

echo "  # # # # # # # # # # # # # # # # # # #
#                                      #
#   Lancement de l'analyse d'exome     #
#                                      #
  # # # # # # # # # # # # # # # # # # #"

# Fonction help
print_usage() {
  printf "How to use :
  -i <path-to-input-files>
  -o <path-to-output>
  -x <path-to-index>
  -g <path-to-gtf>
  "
}

# Valeurs par défaut des paramètres
INPUT="Gautheret/Exome/tp-exome"
OUTPUT="Gautheret/Exome/output"
INDEX="Gautheret/Exome/index"
GTF="Gautheret/gencode.gtf"

# Gestion des paramètres avec des flags
while getopts "i:o:x:g:h" OPTION
do
	case $OPTION in
		i) INPUT=$OPTARG ;;
    o) OUTPUT=$OPTARG ;;
    x) INDEX=$OPTARG ;;
    g) GTF=$OPTARG ;;
    h) print_usage
      exit 1 ;;
	esac
done

# Script

if [ ! -f $OUTPUT ];then
  mkdir $OUTPUT
fi

for file in N T;
# Trimmage des sequences donc elagage des extremites
do
  trimmomatic PE $INPUT/TCRBOA7-$file-WEX-chr16_r1F.fastq.gz $INPUT/TCRBOA7-$file-WEX-chr16_r2F.fastq.gz \
  -baseout $OUTPUT/$file-chr16.fastq  LEADING:20 TRAILING:20 MINLEN:50;

  # Alignement des reads sur le genome de reference
  bwa mem -M -t 4 -A 2 -E 1 $INDEX/chr16.fa.gz \
  $INPUT/TCRBOA7-$file-WEX-chr16_r1F.fastq.gz \
  $INPUT/TCRBOA7-$file-WEX-chr16_r2F.fastq.gz > \
  $OUTPUT/$file-chr16.sam;

  # Conversion de sam vers bam, trie du bam, indexation du bam et calcul de statistiques
  samtools view -b $OUTPUT/$file-chr16.sam -o $OUTPUT/$file-chr16.bam;
  samtools sort $OUTPUT/$file-chr16.bam -o $OUTPUT/$file-chr16-sorted.bam;
  samtools index $OUTPUT/$file-chr16-sorted.bam;
  samtools flagstat $OUTPUT/$file-chr16-sorted.bam > $OUTPUT/$file-chr16-sorted.flagstat;

  #Convert to Mpileup pour avoir un format plus simple pour une comparaison
  samtools mpileup -B -A -f $INDEX/chr16.fa  $OUTPUT/$file-chr16-sorted.bam > $OUTPUT/$file-chr16-mp.mpileup;

done

# Donne les variants entre le tissu normal et le tissu tumoral. Sortie au format VCF
varscan somatic $OUTPUT/N-chr16-mp.mpileup $OUTPUT/T-chr16-mp.mpileup $OUTPUT/chr16.vcf \
  --variants --p-value 0.001 --min-avg-qual 15 --output-vcf 1

# On selectionne les variants somatiques et on convertit la sortie au format BED
grep 'SOMATIC' $OUTPUT/chr16.vcf.snp > $OUTPUT/filtered-chr16.vcf
awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' \
   $OUTPUT/filtered-chr16.vcf > $OUTPUT/chr16.bed

# On reduit le fichier d'annotation GTF aux variants somatiques contenus dans le fichier BED
bedtools intersect -a $GTF -b $OUTPUT/chr16.bed > $OUTPUT/chr16-intersect.txt
# Affichage d'un tableau donnant les genes selectionnes
grep '\sgene\s' $OUTPUT/chr16-intersect.txt | awk '{print " " $1 " " $4 " " $5 " " $16}'

# Produit un fichier d'entree pour Oncotator https://portals.broadinstitute.org/oncotator/
awk '{print $1, "\t", $2, "\t", $2, "\t", $4, "\t",$5}' $OUTPUT/filtered-chr16.vcf > $OUTPUT/chr16.tsv
