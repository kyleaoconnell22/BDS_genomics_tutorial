mkdir resources reference input results input/HG00096 results/HG00096 output_examples
mv input.fastq.zip input/HG0096 | unzip input/HG00096/input.fastq.zip

#Download the indexed reference genome file
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict -o reference/Homo_sapiens_assembly38.dict
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta -o reference/Homo_sapiens_assembly38.fasta
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.alt -o reference/Homo_sapiens_assembly38.fasta.alt
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.amb -o reference/Homo_sapiens_assembly38.fasta.amb
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.ann -o reference/Homo_sapiens_assembly38.fasta.ann
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.bwt -o reference/Homo_sapiens_assembly38.fasta.bwt
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.pac -o reference/Homo_sapiens_assembly38.fasta.pac
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.sa -o reference/Homo_sapiens_assembly38.fasta.sa
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai -o reference/Homo_sapiens_assembly38.fasta.fai

#Download the GATK Resources
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf -o resources/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx -o resources/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz -o resources/1000G_omni2.5.hg38.vcf.gz
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi -o resources/1000G_omni2.5.hg38.vcf.gz.tbi
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz -o resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi -o resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz -o resources/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi -o resources/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf -o resources/Homo_sapiens_assembly38.dbsnp138.vcf
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx -o resources/Homo_sapiens_assembly38.dbsnp138.vcf.idx
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz -o resources/Homo_sapiens_assembly38.known_indels.vcf.gz
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi -o resources/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -o resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi -o resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz -o resources/hapmap_3.3.hg38.vcf.gz
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi -o resources/hapmap_3.3.hg38.vcf.gz.tbi
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/scattered_calling_intervals/ -o resources/scattered_calling_intervals
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list -o resources/wgs_calling_regions.hg38.interval_list

#Download the FASTQ Files from the BDS Bucket
curl https://storage.cloud.google.com/bds-genomics-test-bucket/datatest/HG00096_30xsub/HG00096_CGGACAAC-TCCGGATT_HFHJKDSXX_L004_001.R1.fastq.gz -o input/HG00096/HG00096_CGGACAAC-TCCGGATT_HFHJKDSXX_L004_001.R1.fastq.gz
curl https://storage.cloud.google.com/bds-genomics-test-bucket/datatest/HG00096_30xsub/HG00096_CGGACAAC-TCCGGATT_HFHJKDSXX_L004_001.R2.fastq.gz -o input/HG00096/HG00096_CGGACAAC-TCCGGATT_HFHJKDSXX_L004_001.R2.fastq.gz

#Download example outputs from BDS Bucket
curl https://storage.cloud.google.com/bds-genomics-test-bucket/output_examples/HG00096_deepvariant.visual_report_example.html -o output_examples/HG00096_deepvariant.visual_report_example.html
curl https://storage.cloud.google.com/bds-genomics-test-bucket/output_examples/snps.recal -o output_examples/snps.recal
curl https://storage.cloud.google.com/bds-genomics-test-bucket/output_examples/snps.recal.idx -o output_examples/snps.recal.idx
curl https://storage.cloud.google.com/bds-genomics-test-bucket/output_examples/snps.tranches -o output_examples/snps.tranches
curl https://storage.cloud.google.com/bds-genomics-test-bucket/output_examples/HG00096_haplotypecaller.vcf -o output_examples/HG00096_haplotypecaller.vcf 

