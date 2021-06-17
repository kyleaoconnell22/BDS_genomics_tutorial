#!/usr/bin/python
import sys, argparse, subprocess, time, os, re, glob, numpy as np
from ftplib import FTP

# Globals
path_to_data = ''
path_to_ref = ''
path_to_results = os.path.dirname(os.path.abspath(__file__))
index_genome = False
verbosity = False
threads = 1
deleteResultFiles = False


def create_directory(mypath, verbose=False):
    try:
        if not os.path.exists(mypath):
            os.makedirs(mypath)
            if verbose:
                print("Directory created: ", mypath)
        else:
            if verbose:
                print("Directory already exists: ", mypath)

    except OSError as exc:  # Guard against race condition
        if exc.errno != errno.EEXIST:
            print("Directory could not be created: ", path)
            raise


def run_command(cmd, verbose=False, dry_run=False):
    if verbose:
        print(f'Executing command: {cmd}...', end='', flush=True)

    if dry_run:
        if verbose:
            print('Dry run: nothing executed.', end='')
    else:
        # https://docs.python.org/3/library/subprocess.html
        if type(cmd) == str:
            result = subprocess.run(cmd.split(' '), stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
        elif type(cmd) == list:
            result = subprocess.run(cmd, stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
        else:
            print(f"Invalid type: {type(cmd)}")

    if verbose:
        print('done.', flush=True)
        print("Result is: ", result)
    return result


def run_command_out_to_file(cmd, outfile, verbose=False, dry_run=False):
    if verbose:
        print(f'Executing command: {cmd} > {outfile}', end='', flush=True)

    if dry_run:
        if verbose:
            print('Dry run: nothing executed.', end='')
    else:
        # https://docs.python.org/3/library/subprocess.html
        f = open(outfile, "w")
        if type(cmd) == str:
            result = subprocess.run(cmd.split(' '), stdout=f)
        elif type(cmd) == list:
            result = subprocess.run(cmd, stdout=f)
        else:
            print(f"Invalid type: {type(cmd)}")
    if verbose:
        print('done.', flush=True)
        print("Result is: ", result)
    return result


def run_2pipedcommands(cmd, cmd2, verbose=False, dry_run=False):
    if verbose:
        print(f'Executing command: {cmd} | {cmd2} ...', end='', flush=True)

    if dry_run:
        if verbose:
            print('Dry run: nothing executed.', end='')
    else:
        # result = subprocess.run(cmd.split(' '), stdout=subprocess.PIPE)
        p1 = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cmd2.split(' '), stdin=p1.stdout, stdout=subprocess.PIPE)
        result = p2.communicate()[0]

    if verbose:
        print('done.', flush=True)
        print("Result is: ", result)
    return result


def find_mean(file):
    d = []
    fh = open(file, 'r')
    for line in fh:
        line = line.strip()
        d.append(int(line.split('\t')[2]))
    result = round(np.mean(d),3)
    print ('Mean coverage = ', result)
    return result


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("path_to_data", help="Absolute path to FASTQ data (data/subject/files)")
    parser.add_argument("path_to_ref", help="Absolute path for reference genome")
    parser.add_argument("--path_to_results", help="Absolute path for results (default is current directory)")
    parser.add_argument("--path_to_resources", help="Absolute path for GATK resources dir)")
    parser.add_argument("--variant_caller", help="Variant caller can be either haplotypecaller or deepvariant")
    parser.add_argument("--filter_vcf", help="VQSR, CNN, or NONE")
    parser.add_argument("--index_genome", help="Index the genome? (True/False)", action='store_true')
    parser.add_argument("--verbose", help="Print debugging information", action='store_true')
    parser.add_argument("--threads", help="Number of threads")
    parser.add_argument("--deleteResultFiles", help="Delete results directory after execution", action='store_true')

    args = parser.parse_args()
    if args.path_to_data:
        path_to_data = args.path_to_data
    if args.path_to_ref:
        path_to_ref = args.path_to_ref
    if args.path_to_resources:
        path_to_resources=args.path_to_resources
    if args.path_to_results:
        path_to_results = args.path_to_results
    if args.index_genome:
        index_genome = args.index_genome
    if args.verbose:
        verbosity = args.verbose
    if args.threads:
        threads = int(args.threads)
    if args.deleteResultFiles:
        deleteResultFiles = args.deleteResultFiles
    if args.variant_caller:
        var_caller = args.variant_caller
    if args.filter_vcf:
        filter=args.filter_vcf

    print('Data path is: ', path_to_data)
    print('Ref  path  is: ', path_to_ref)
    print('Resources path is: ', path_to_resources)
    print('Index Genome option is: ', index_genome)
    print('Variant Caller is: ', var_caller)
    print('Verbosity is: ', verbosity)
    print("Number of threads is: ", threads)
    print("Filter VCF option is: ",filter)
    print("deleteResultFiles option is: ", deleteResultFiles)

    if index_genome:
        run_command(f'bwa index {path_to_ref}', verbose=verbosity)

    # Traverse FASTQ data and look for paired reads only
    start_time = time.time()
    subjects = os.listdir(path_to_data)

    # If any subjects create results directory if it doesn't exist
    if subjects:
        create_directory(path_to_results, verbose=verbosity)

    for s in subjects:
        print('sample=',s)
        fastq_files = os.listdir(os.path.join(path_to_data, s))
        #print('all_files: ',fastq_files)
        # Look for _1 FASTQ files
        p = []
        for fastq in fastq_files:
                if fastq.endswith('.R1.fastq.gz'):
                        p.append(fastq)
        files = [os.path.join(path_to_data, s, f) for f in p]
        i=1
        # if we have files create an appropriate subdirectory
        if files:
            create_directory(os.path.join(path_to_results, s), verbose=verbosity)

        # For each paired read
        for f1 in files:
            # Confirm _2 FASTQ file exists
            f2 = f1.replace("R1", "R2")
            if not os.path.exists(f2):
                print("Matching Paired read file missing: ", f2)
                continue

            # If it does, run mapping commands
            finalBAM = f1.replace("R1", "_12")
            finalBAM = finalBAM.replace(path_to_data, path_to_results)
            if "fastq.gz" in finalBAM:
                finalBAM = finalBAM.replace(".fastq.gz", ".bam")
            else:
                finalBAM = finalBAM.replace(".fastq", ".bam")

            # Create our BAM
            #read group info
            LBstr='lib'+str(i)
            #LB=library# which should be unique for each because they are PCR free
            seq = f1.split('/')[-1]
            base = os.path.basename(f1)
            IDstr=base.split('.')[0].split('_')[2]+'.'+seq.split('.')[0].split('_')[3] #flowcell.lane
            PUstr=IDstr+'.'+base.split('.')[0].split('_')[1].replace('-','_') #flowcell.lane.barcode
            rg_string = f'@RG\\tID:{IDstr}\\tLB:{LBstr}\\tPL:LLUMINA\\tPU:{PUstr}\\tSM:{s}'
            print('read group string=',rg_string)
            bwamem = ['bwa', 'mem', '-Y', '-K', '10000000', '-t', str(threads), '-R', rg_string,
                      path_to_ref, f1, f2]
            run_command_out_to_file(bwamem, "out.sam", verbose=verbosity)

            # Sort Sam File
            run_command(f'gatk SortSam -I out.sam -O {finalBAM} -SO coordinate', verbose=verbosity)
            i = i+1

        # Combine all BAMs for this subject 's'
        combinedBAM = os.path.join(path_to_results, s, s + '.bam')
        combinedVCF = os.path.join(path_to_results, s, s + '_' + var_caller + '.vcf')
        filtVCF = os.path.join(path_to_results, s, s + '_' + var_caller + '.filt.vcf')
        allBAMs = glob.glob(os.path.join(path_to_results, s, '*.bam'))
        BAMarguments = ' '.join([f'I={x}' for x in allBAMs])
        dupBAM = os.path.join(path_to_results, s, s + '.markdup.bam')
        indexBAM = os.path.join(path_to_results, s, s + '.markdup.bai')
        GIAB_VCF = os.path.join(path_to_resources, 'NA12878-NISTv2.19v2_sorted.vcf')
        GIAB_BED = os.path.join(path_to_resources, 'NA12878-NISTv2.19.sorted.bed')
        known_sites = os.path.join(path_to_resources, 'Homo_sapiens_assembly38.known_indels.vcf.gz')
        training_data1 = os.path.join(path_to_resources, '1000G_omni2.5.hg38.vcf.gz')
        training_data2 = os.path.join(path_to_resources, 'hapmap_3.3.hg38.vcf.gz')

        # Mark Duplicates - use MULTIPLE INPUTS if multiple files present
        run_command(f'gatk MarkDuplicates {BAMarguments} O={dupBAM} M=metrics.txt',
                    verbose=verbosity)

        if var_caller == 'haplotypecaller':
            # BQSR Report
            run_command(
                f'gatk BaseRecalibrator -I {dupBAM} -O bqsr_out.txt --known-sites {known_sites} -R {path_to_ref}',
                verbose=verbosity)

            # Rescale bam file
            run_command(
                f'gatk ApplyBQSR -R {path_to_ref} -I {dupBAM} -O {combinedBAM} -bqsr-recal-file bqsr_out.txt',verbose=verbosity)

            # Run haplotype caller
            run_command_out_to_file(f'gatk HaplotypeCaller -I {combinedBAM} -O {combinedVCF} -R {path_to_ref} --native-pair-hmm-threads {str(threads)} --output-mode EMIT_VARIANTS_ONLY','out.vcf',verbose=verbosity)

            # See how many reads were successfully mapped
            run_command(f'samtools flagstat {combinedBAM}', verbose=verbosity)

            # Calculate mapping read depth
            run_command_out_to_file(f'samtools depth {combinedBAM}', "depth.out", verbose=verbosity)
            find_mean("depth.out")

            # Delete output files
            if deleteResultFiles:
                run_command(f'rm -rf {path_to_results}')
            s
            if filter == 'VQSR':
                # Run VQSR SNPs
                run_command(f'gatk VariantRecalibrator -V {combinedVCF} --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 -an QD -an FS -an MQ -an SOR -an DP -mode SNP -resource:hapmap,known=false,training=true,truth=true,prior=15 resources/hapmap_3.3.hg38.vcf.gz -resource:omni,known=false,training=true,truth=true,prior=12 resources/1000G_omni2.5.hg38.vcf.gz -resource:1000G,known=false,training=true,truth=false,prior=10 resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz -resource:dbsnp,known=true,training=false,truth=false,prior=7 resources/Homo_sapiens_assembly38.dbsnp138.vcf -O snps.recal --tranches-file snps.tranches',verbose=verbosity)

                # Run VQSR indels
                run_command(f'gatk VariantRecalibrator -V {combinedVCF} -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -mode INDEL -resource:mills,known=false,training=true,truth=true,prior=12 resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O indels.recal --tranches-file indels.tranches -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP',verbose=verbosity)

                # Apply VQSR
                run_command(f'gatk ApplyVQSR -R {path_to_ref} -V {combinedVCF} -O {filtVCF} --tranches-file snps.tranches --recal-file snps.recal -mode SNP', verbose=verbosity)

                #print ('Number of SNPs in final vcf file after VQSR:  ')
                run_command(f'vcftools --vcf {combinedVCF}', verbose=verbosity)

            elif filter == 'NONE':
                pass

            elif filter == 'CNN':
                print('This option is not ready yet, please use VQSR or NONE')
                    #can also add gatk CNN filter model but it needs some more troubleshooting
                #https://github.com/gatk-workflows/gatk4-cnn-variant-filter
                #https://gatk.broadinstitute.org/hc/en-us/articles/360035530952

        # Run Deepvariant
        elif var_caller == 'deepvariant':
            # Index BAM file
            print ('Indexing BAM')
            run_command(f'gatk BuildBamIndex -I {dupBAM} -O {indexBAM}',verbose=verbosity)

            # Run Deepvariant
            print ('Running Deepvariant')
            run_command('sh run-deepvariant.sh', verbose=verbosity)

            #print ('Number of SNPs in final vcf file after DeepVariant:  ')
            run_command(f'vcftools --vcf {combinedVCF}', verbose=verbosity)

            #Summarize BAM file
            # See how many reads were successfully mapped
            run_command(f'samtools flagstat {dupBAM}', verbose=True)

            # Calculate mapping read depth
            run_command_out_to_file(f'samtools depth {dupBAM}', "depth.out", verbose=verbosity)

            # Run function to find average of depth/read across all reads
            find_mean("depth.out")

            # Delete output files
            if deleteResultFiles:
                run_command(f'rm -rf {path_to_results}')

        else:
            print('variant_caller command accepts either haplotypecaller or deepvariant as options')
            pass

        #compare to GIAB
        #if s == 'HG001':
        #    run_command('sh run-hap.py.sh', verbose=verbosity)

elapsed_time = time.time() - start_time
print(f"Time taken was {round(elapsed_time)} seconds.")
