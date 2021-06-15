mkdir resources reference input results input/HG00096 results/HG00096 output_examples
unzip input.fastq.zip
mv *.gz input/HG00096/

#Download the indexed reference genome file
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict -o reference/Homo_sapiens_assembly38.dict
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta -o reference/Homo_sapiens_assembly38.fasta
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt -o reference/Homo_sapiens_assembly38.fasta.64.alt
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb -o reference/Homo_sapiens_assembly38.fasta.64.amb
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann -o reference/Homo_sapiens_assembly38.fasta.64.ann
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt -o reference/Homo_sapiens_assembly38.fasta.64.bwt
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac -o reference/Homo_sapiens_assembly38.fasta.64.pac
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa -o reference/Homo_sapiens_assembly38.fasta.64.sa
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai -o reference/Homo_sapiens_assembly38.fasta.fai

#Download the GATK Resources
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz -o resources/Homo_sapiens_assembly38.known_indels.vcf.gz
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi -o resources/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi

