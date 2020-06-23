### Use bowtie2 to align reads

```
bowtie2 bowtie2-indexes/dm6 -p 8 --end-to-end --very-sensitive --no-mixed --no-discordant -q --phred33 -I 10 -X 700 \
    -1 R1.fastq30 R2.fastq30 | samtools view -bT dm6.fa - | samtools sort - aln.bam
```

### Index the bam file

```
samtools index aln.bam
```

