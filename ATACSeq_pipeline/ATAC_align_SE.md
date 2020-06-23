### Use bowtie1 to align reads

```
bowtie bowtie-indexes/dm6 -p 8 -v 2 -k 1 -m 3 -t --sam-nh --best -y --strata -q --sam \
    R.fastq30 | samtools view -bT dm6.fa - | samtools sort - aln.bam
```

### Index the bam file

```
samtools index aln.bam
```
