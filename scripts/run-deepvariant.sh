BIN_VERSION="1.1.0"  
S='HG00096'  
INPUT_DIR="${PWD}/results/${S}"  
REF_DIR="${PWD}/reference"  
REF="Homo_sapiens_assembly38.fasta"  
BAM="${S}.markdup.bam"  
OUTPUT_DIR="${INPUT_DIR}"  
OUTPUT_VCF="${S}_deepvariant.vcf"  

sudo docker run -v "${INPUT_DIR}":"/input" -v "${OUTPUT_DIR}":"/output" -v "${REF_DIR}":"/refdir" google/deepvariant:"1.1.0" /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref="/refdir/${REF}"  --reads="/input/${BAM}"  --output_vcf="/output/${OUTPUT_VCF}" --num_shards=4 --intermediate_results_dir /output/intermediate_results_dir  
