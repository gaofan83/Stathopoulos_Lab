### 1. HOMER peak calling in sensitive mode

#### 1.1 make HOMER tags
```
homer-4.7/bin/makeTagDirectory SampleHomerTags Sample_nochrM.bam
```

#### 1.2 call peaks
```
homer-4.7/bin/findPeaks SampleHomerTags -localSize 50000 -minDist 50 -size 150 -fragLength 0 -o SamplelS50000mD50s150fL0 2> HOMER.err
grep 000 SamplelS50000mD50s150fL0 | grep chr - | grep -v = | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"225"\t"$5}' - | sort -k 1d,1 -k 2n,2 > SamplelS50000mD50s150fL0.whole.bed
```

#### 1.3 remove peaks overlapping blacklisted regions
```
intersectBed -a SamplelS50000mD50s150fL0.whole.bed -b blacklist.bed -v > SamplelS50000mD50s150fL0.bed
```

### 2. model the background

#### 2.1 generate coverage signal (RPM) file
```
deepTools-2.4.2_develop/bin/bamCoverage -b Sample_nochrM.bam -of bedgraph -bs 1 --scaleFactor 1 --skipNonCoveredRegions --smoothLength 50 -o Sample_nochrM.count.bg4
```

#### 2.2 extend peaks and extract background regions.
```
bedtools random -l 150 -seed 1 -n 100000 -g chrom.sizes | bedtools intersect -a - -b <(bedtools slop -b 9925 -i SamplelS50000mD50s150fL0.whole.bed -g chrom.sizes) -v | bedtools sort -i - > SamplelS50000mD50s150fL0.background.bed
```

#### 2.3 calculate mean and std of the background noise
```
SamplelS50000mD50s150fL0backgroundscores=$(bedtools intersect -u -a Sample_nochrM.count.bg4 -b SamplelS50000mD50s150fL0.background.bed | awk '{if($4!=0){array1+=1;arrayX+=log($4)/log(2);arrayXsq+=(log($4)/log(2))^2}} END {print arrayX/array1" "sqrt(arrayXsq/array1 - arrayX^2/array1^2)}')
SamplelS50000mD50s150fL0backgroundmean=$(echo $SamplelS50000mD50s150fL0backgroundscores | cut -d ' ' -f1)
SamplelS50000mD50s150fL0backgroundstd=$(echo $SamplelS50000mD50s150fL0backgroundscores | cut -d ' ' -f2)
```

### 3. normalize ATAC signal
```
awk '{print $1"\t"$2"\t"$3"\t"2^((log($4)/log(2)-'$SamplelS50000mD50s150fL0backgroundmean')/'$SamplelS50000mD50s150fL0backgroundstd')-1}' Sample_nochrM.count.bg4 > Sample_nochrM.normalized.bg4
wigToBigWig -clip Sample_nochrM.normalized.bg4 dm3.chrom.sizes Sample_nochrM.normalized.bigWig
```
