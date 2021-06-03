#!/bin/bash

#Get the experiment name and output file name
status_code=$(curl --write-out %{http_code} --silent --output /dev/null http://metadata.google.internal/computeMetadata/v1/instance/attributes/experimentname -H "Metadata-Flavor: Google")
echo "Status code for threads tag check was $status_code"
if [[ "$status_code" == 200 ]] ; then
EXPERIMENTNAME=$(curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/experimentname -H "Metadata-Flavor: Google")
else
EXPERIMENTNAME=$(date -u +%y_%m_%d_%H%M%S)
fi
FILENAME="output_$EXPERIMENTNAME.log"
echo "Experiment Name is $EXPERIMENTNAME"
echo "Output filename is $FILENAME"

#Get the threads
status_code=$(curl --write-out %{http_code} --silent --output /dev/null http://metadata.google.internal/computeMetadata/v1/instance/attributes/mappingthreads -H "Metadata-Flavor: Google")
echo "Status code for threads tag check was $status_code"
if [[ "$status_code" == 200 ]] ; then
  THREADS=$(curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/mappingthreads -H "Metadata-Flavor: Google") 
else
  THREADS=8
fi
echo "Number of threads to use is: $THREADS"

#Get the variant caller
status_code=$(curl --write-out %{http_code} --silent --output /dev/null http://metadata.google.internal/computeMetadata/v1/instance/attributes/variantcaller -H "Metadata-Flavor: Google")
echo "Status code for threads tag check was $status_code"
if [[ "$status_code" == 200 ]] ; then
  VAR_CALLER=$(curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/variantcaller -H "Metadata-Flavor: Google") 
else
  VAR_CALLER=haplotypecaller
fi
echo "Variant Caller is : $VAR_CALLER"

#Get the variant caller
status_code=$(curl --write-out %{http_code} --silent --output /dev/null http://metadata.google.internal/computeMetadata/v1/instance/attributes/filter -H "Metadata-Flavor: Google")
echo "Status code for threads tag check was $status_code"
if [[ "$status_code" == 200 ]] ; then
  filter=$(curl http://metadata.google.internal/computeMetadata/v1/instance/attributes/filter -H "Metadata-Flavor: Google") 
else
  filter=VQSR
fi
echo "Filter options is : $filter"

echo "Copying scripts from storage bucket..."
cd /home/usa_kyoconnell_deloitte_com/
gsutil cp gs://bds-genomics-test-bucket/scripts/execute_mapping_v2.4.2.py .
gsutil cp  gs://bds-genomics-test-bucket/scripts/run-deepvariant.sh .
echo "Copy complete"
echo 'home dir' $HOME 

#init conda
export PATH="$HOME/miniconda3/bin:$PATH"
cp /home/usa_kyoconnell_deloitte_com/miniconda3/bin/vcftools .
echo "Starting mapping job"
#/usr/bin/python3 execute_mapping_v2.4.1.py /home/usa_kyoconnell_deloitte_com/datatest_sub30x /home/usa_clobb_deloitte_com/BioProject/refsBWAKIT/hs38DH.fa --path_to_results /home/usa_clobb_deloitte_com/BioProject/results --threads $THREADS --variant_caller $VAR_CALLER > $FILENAME 2>&1
#/usr/bin/python3 execute_mapping_v2.4.1.py /home/usa_clobb_deloitte_com/BioProject/datatest_30xNYCOneSubj /home/usa_clobb_deloitte_com/BioProject/refsBWAKIT/hs38DH.fa --path_to_results /home/usa_clobb_deloitte_com/BioProject/results --threads $THREADS --variant_caller $VAR_CALLER > $FILENAME 2>&1
/usr/bin/python3 execute_mapping_v2.4.2.py /home/usa_kyoconnell_deloitte_com/datatest_sub30x /home/usa_kyoconnell_deloitte_com/reference/Homo_sapiens_assembly38.fasta  --path_to_results /home/usa_kyoconnell_deloitte_com/results --threads $THREADS --variant_caller $VAR_CALLER --filter_vcf $filter --path_to_resources /home/usa_kyoconnell_deloitte_com/resources > $FILENAME 2>&1

echo "Job complete"

#Upload to gs
gsutil cp $FILENAME gs://bds-genomics-test-bucket/output_logs/
gsutil cp /home/usa_kyoconnell_deloitte_com/results/*/*.vcf gs://bds-genomics-test-bucket/outputVCFs/
gsutil cp /home/usa_kyoconnell_deloitte_com/results/*/*.html gs://bds-genomics-test-bucket/outputVCFs/

#Terminate
sudo shutdown -h +1
