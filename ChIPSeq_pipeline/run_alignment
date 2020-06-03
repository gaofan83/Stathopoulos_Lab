chr_mouse="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX"
chr_human="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chrX"
chr_drosophila="chr2L chr2R chr3L chr3R chr4 chrX"
#concatenate files 
while read line;
do
java -jar $TRIMMOMATIC_DIR/trimmomatic-0.38.jar SE -threads 32 -phred33 $line.combined.fastq.gz $line.combined.trim.fastq.gz \
      ILLUMINACLIP:$TRIMMOMATIC_DIR/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36

bowtie2 -x /home/fgao/reference_genome/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome -U $line.combined.trim.fastq.gz -S $line.sam -p 32
samtools view -bS $line.sam > $line.bam
rm $line.sam
samtools sort $line.bam -o $line.sort.bam
samtools rmdup -s $line.sort.bam $line.rmdup.bam
samtools idxstats $line.sort.bam > $line.mapping.txt
samtools idxstats $line.rmdup.bam > $line.mapping_rmdup.txt
done < sample_sel.txt

while read line;
do
samtools index $line.sort.bam
samtools index $line.rmdup.bam
samtools view -b $line.rmdup.bam $chr_drosophila > $line.rmdup.nochrM.bam
samtools index $line.rmdup.nochrM.bam
samtools idxstats $line.sort.bam | awk '{s+=$3+$4} END {print "Total\t",s}' > $line.mapping.txt
samtools idxstats $line.sort.bam | awk '{s+=$3} END {print "Mapped\t",s}' >> $line.mapping.txt
samtools idxstats $line.sort.bam | awk '{if($1=="chrM") print "chrM\t",$3}' >> $line.mapping.txt
samtools idxstats $line.rmdup.bam | awk '{s+=$3} END {print "Mapped_rmdup\t",s}' >> $line.mapping.txt
samtools idxstats $line.rmdup.nochrM.bam | awk '{s+=$3} END {print "Mapped_rmdup_nochrM\t",s}' >> $line.mapping.txt
done < sample_ID.txt
