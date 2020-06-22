### 1. Check reads quality using fastqc

```
fastqc R1.fastq -o FastQCk6R1 -k 6 && fastqc R2.fastq -o FastQCk6R2 -k 6
```

### 2. Use cutadapt to trim common adaptors at 5' and 3' ends

```
cutadapt -a CTGTCTCTTATACAC R1.fastq | cutadapt -a CGTATGCCGTCTTCTGCTTG - | cutadapt -g TGCCGTCTTCTGCTTG - | cutadapt \
    -g GGTAACTTTGTGTTT - | cutadapt -g CTTTGTGTTTGA - | cutadapt -a CACTCGTCGGCAGCGTTAGATGTGTATAAG - | cutadapt \
    -a GAAGAGCACACGTCTGAACTCC - > R1.trim1.fastq && \
cutadapt -a CTGTCTCTTATACAC R2.fastq | cutadapt -a CGTATGCCGTCTTCTGCTTG - | cutadapt -g TGCCGTCTTCTGCTTG - | cutadapt \
    -g GGTAACTTTGTGTTT - | cutadapt -g CTTTGTGTTTGA - | cutadapt -a CACTCGTCGGCAGCGTTAGATGTGTATAAG - | cutadapt \
    -a GAAGAGCACACGTCTGAACTCC - > R2.trim1.fastq
```

### 3. Use Trimmomatic to trim indexed barcodes
The adaptor sequence file can be found [here](https://github.com/brianpenghe/Stathopoulos_Lab/blob/master/ATACSeq_pipeline/AllAdaptors.fa)

```
trimmomatic-0.33.jar PE -threads 4 R1.trim1.fastq R2.trim1.fastq R1.trim2.fastq R1.trim2.unpaired.fastq \
    R2.trim2.fastq R2.trim2.unpaired.fastq \
    ILLUMINACLIP:Alladaptors.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
```

### 4. Check reads quality again

```
fastqc R1.trim2.fastq -o FastQCk6R1 -k 6 && fastqc R2.trim2.fastq -o FastQCk6R2 -k 6
```

### 5. Trim reads down to 30bp
The script can be found [here](https://github.com/brianpenghe/fastq-scripts/blob/master/trimfastq.py)
```
python2 trimfastq.py R1.trim2.fastq 30 -stdout > R1.fastq30 && \
python2 trimfastq.py R2.trim2.fastq 30 -stdout > R2.fastq30
```
