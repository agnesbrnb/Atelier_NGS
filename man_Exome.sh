#bwa index -a bwtsw chr16.fa.gz

# Initialisation et erreurs


echo "  # # # # # # # # # # # # # # # # # # #
#                                      #
#   Lancement de l'analyse d'exome     #
#                                      #
  # # # # # # # # # # # # # # # # # # #"

# Script

INPUT=${1:-"Gautheret/Exome/tp-exome"}
OUTPUT=${2:-"Gautheret/Exome"}
INDEX=${3:-"Gautheret/Exome/index"}

if [ ! -f $OUTPUT/bwa ];then
  mkdir $OUTPUT/bwa
fi

for file in N T;
do trimmomatic PE $INPUT/TCRBOA7-$file-WEX-chr16_r1F.fastq.gz $INPUT/TCRBOA7-$file-WEX-chr16_r2F.fastq.gz \
-baseout $OUTPUT/$file-chr16.fastq  LEADING:20 TRAILING:20 MINLEN:50;

bwa mem -M -t 4 -A 2 -E 1 $INDEX/chr16.fa.gz \
$INPUT/TCRBOA7-$file-WEX-chr16_r1F.fastq.gz \
$INPUT/TCRBOA7-$file-WEX-chr16_r2F.fastq.gz > \
$OUTPUT/bwa/$file-chr16.sam;

samtools view -b $OUTPUT/bwa/$file-chr16.sam -o $OUTPUT/bwa/$file-chr16.bam;
samtools sort $OUTPUT/bwa/$file-chr16.bam -o $OUTPUT/bwa/$file-chr16-sorted.bam;
samtools index $OUTPUT/bwa/$file-chr16-sorted.bam;
samtools flagstat $OUTPUT/bwa/$file-chr16-sorted.bam > $OUTPUT/bwa/$file-chr16-sorted.flagstat;

#Convert to Mpileup
samtools mpileup -B -A -f $INDEX/chr16.fa  $OUTPUT/bwa/$file-chr16-sorted.bam > $OUTPUT/bwa/$file-chr16-mp.mpileup;

done
