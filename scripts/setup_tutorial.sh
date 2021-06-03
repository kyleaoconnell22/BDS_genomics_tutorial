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
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz -o resources/Homo_sapiens_assembly38.known_indels.vcf.gz
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi -o resources/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi

