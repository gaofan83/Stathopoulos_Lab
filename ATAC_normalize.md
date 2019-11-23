#HOMER
homer-4.7/bin/makeTagDirectory SampleHomerTags Sample_nochrM.bam
homer-4.7/bin/findPeaks SampleHomerTags -localSize 50000 -minDist 50 -size 150 -fragLength 0 -o SamplelS50000mD50s150fL0 2> HOMER.err
grep 000 SamplelS50000mD50s150fL0 | grep chr - | grep -v = | awk '{print $2"\t"$3"\t"$4"\t"$1"\t"225"\t"$5}' - | sort -k 1d,1 -k 2n,2 > SamplelS50000mD50s150fL0.whole.bed
intersectBed -a SamplelS50000mD50s150fL0.whole.bed -b blacklist.bed -v > SamplelS50000mD50s150fL0.bed
deepTools-2.4.2_develop/bin/bamCoverage -b Sample_nochrM.bam -of bedgraph -bs 1 --scaleFactor 1 --skipNonCoveredRegions --smoothLength 50 -o Sample_nochrM.count.bg4
bedtools random -l 150 -seed 1 -n 100000 -g chrom.sizes | bedtools intersect -a - -b <(bedtools slop -b 9925 -i SamplelS50000mD50s150fL0.whole.bed -g chrom.sizes) -v | bedtools sort -i - > SamplelS50000mD50s150fL0.background.bed
SamplelS50000mD50s150fL0backgroundscores=$(bedtools intersect -u -a Sample_nochrM.count.bg4 -b SamplelS50000mD50s150fL0.background.bed | awk '{if($4!=0){array1+=1;arrayX+=log($4)/log(2);arrayXsq+=(log($4)/log(2))^2}} END {print arrayX/array1" "sqrt(arrayXsq/array1 - arrayX^2/array1^2)}')
SamplelS50000mD50s150fL0backgroundmean=$(echo $SamplelS50000mD50s150fL0backgroundscores | cut -d ' ' -f1)
SamplelS50000mD50s150fL0backgroundstd=$(echo $SamplelS50000mD50s150fL0backgroundscores | cut -d ' ' -f2)
awk '{print $1"\t"$2"\t"$3"\t"2^((log($4)/log(2)-'$SamplelS50000mD50s150fL0backgroundmean')/'$SamplelS50000mD50s150fL0backgroundstd')-1}' Sample_nochrM.bamp.count.bg4 > Sample_nochrM.normalized.bg4
