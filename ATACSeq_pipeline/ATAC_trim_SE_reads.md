### 1. Check reads quality using fastqc

```
fastqc R.fastq -o FastQCk6
```

### 2. Use cutadapt to trim common adaptors at 5' and 3' ends

```
cutadapt -a CTGTCTCTTATACAC R.fastq | cutadapt -a CGTATGCCGTCTTCTGCTTG - | cutadapt -g TGCCGTCTTCTGCTTG - | cutadapt \
    -g GGTAACTTTGTGTTT - | cutadapt -g CTTTGTGTTTGA - | cutadapt -a CACTCGTCGGCAGCGTTAGATGTGTATAAG - | cutadapt \
    -a GAAGAGCACACGTCTGAACTCC - > R.trim1.fastq 
```

### 3. Use Trimmomatic to trim indexed barcodes
The adaptor sequence file can be found [here](https://github.com/brianpenghe/Stathopoulos_Lab/blob/master/ATACSeq_pipeline/AllAdaptors.fa)

```
trimmomatic-0.33.jar SE -threads 4 R.trim1.fastq R.trim2.fastq \
    ILLUMINACLIP:Alladaptors.fa:2:30:10 MAXINFO:35:0.9 MINLEN:30
```

### 4. Check reads quality again

```
fastqc R.trim2.fastq -o FastQCk6R -k 6
```

### 5. Trim reads down to 30bp
The script can be found [here](https://github.com/brianpenghe/fastq-scripts/blob/master/trimfastq.py)
```
python2 trimfastq.py R.trim2.fastq 30 -stdout > R.fastq30
```
