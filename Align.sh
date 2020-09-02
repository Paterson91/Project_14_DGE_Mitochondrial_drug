for i in `find . -name "*_R1_sorted.fastq.gz" | sed 's/_R1_sorted.fastq.gz//'`
	do
		echo ""
		echo ">>>>>>>>>>>>>>>>>>>> $i 1st Pass Alignment begun <<<<<<<<<<<<<<<<<<<<"
		echo ""
		STAR --genomeDir /home/ap14958/STAR_Indexes/Homo_sapiens_GRCh38/ \
    --runMode alignReads \
		--readFilesIn $i\_R1_sorted.fastq.gz $i\_R2_sorted.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix STAR_output/`basename $i` \
    --outSJfilterReads Unique \
		--outSAMtype BAM SortedByCoordinate --runThreadN 100
		echo ""
		echo ">>>>>>>>>>>>>>>>>>>> $i 1st Pass Alignment complete <<<<<<<<<<<<<<<<<<<<"
		echo ""
	done
