featureCounts \
-p -g gene_name -t exon \
-a /home/ap14958/GTF/Homo_sapiens/Homo_sapiens.GRCh38.98.gtf \
-o Counts_stranded_gtf.gene -F GTF -T 64 \
-s 0 \
*.bam

tail -n +2 Counts_stranded_gtf.gene | cut -f1,7- | sed 's/Aligned.sortedByCoord.out.bam//g'> Counts_stranded_gtf.gene_tidied
