source ~/.bashrc

nthreads=1

outdir=output/S1_$(date +%C%y%m%d_%H_%M_%S)/
indir=raw_data/2017-07-05_GCMP-EMP-shipment1/
sample2bc=raw_data/EMP_GMCP_samples_sequenced_20170619.txt

mkdir -p ${outdir}

awk -F '\t' -v outdir=${outdir} 'NR>1 && $5=="GCMP1_Lane7_S7_L007" {print $1"\t"$2"\t""GCMP1" > outdir"/bcalls_GCMP1.txt"} NR>1 && $5=="GCMP2_Lane8_S8_L008" {print $1"\t"$2"\tGCMP2" > outdir"/bcalls_GCMP2.txt"}' ${sample2bc}


mkdir -p ${outdir}/demultiplexed/bcalls_split/
mkdir -p ${outdir}/demultiplexed/fastqs

fwd=(${indir}/GCMP1_Lane7_S7_L007_R1_001.fastq.gz ${indir}/GCMP2_Lane8_S8_L008_R1_001.fastq.gz)
rev=(${indir}/GCMP1_Lane7_S7_L007_R2_001.fastq.gz ${indir}/GCMP2_Lane8_S8_L008_R2_001.fastq.gz)
ind=(${indir}/GCMP1_Lane7_S7_L007_I1_001.fastq.gz ${indir}/GCMP2_Lane8_S8_L008_I1_001.fastq.gz)
bcalls=(${outdir}/bcalls_GCMP1.txt ${outdir}/bcalls_GCMP2.txt)

for run in 0 1; do

	split -l 100 ${bcalls[run]} ${outdir}/demultiplexed/bcalls_split/bcalls_${run}_

	for file in ${outdir}/demultiplexed/bcalls_split/bcalls_${run}_*; do

		fastq-multx -B ${file} -m 1.5 -d 0 ${ind[run]} ${fwd[run]} ${rev[run]} -o n/a -o ${outdir}/demultiplexed/fastqs/%_R1.fastq.gz -o ${outdir}/demultiplexed/fastqs/%_R2.fastq.gz

		rm ${outdir}/demultiplexed/fastqs/unmatched_*.fastq.gz

	done

done


mkdir -p ${outdir}/filtered/
mkdir -p ${outdir}/shortreads/merged/

for file in ${outdir}/demultiplexed/fastqs/*_R1.fastq.gz; do

	# get Sample IDs that match mapping file
	filename=$(basename $file)
	sampleid=${filename/_R1.fastq.gz/}

	## if the file is not empty:
	if [ -s ${file} ]

	then

        ## concatenate the forward read to the reverse complement of the reverse read. might consider trimming the reverse read first if they are low quality?
		paste <(zcat ${file}) <(vsearch --fastx_revcomp ${file/R1/R2} --fastqout -) |

            ## fix headers so they contain the sequence itself and the sample ID
            awk '(NR % 4) == 1 {print ">"} (NR % 4) == 2 {print $1$2}' |
                awk -v sample=${sampleid} '$0 !~ ">" && length($0) > 0 {print ">"$0";sample="sample"\n"$0}' |

                    ## derep sequences
                    vsearch --derep_fulllength - --sizeout --output - |

                        ## denoise sequences
                        vsearch --cluster_unoise - --sizein --sizeout --fasta_width 0 --centroids ${outdir}/filtered/${sampleid}_denoised.fasta

        ## remove chimeras
        vsearch --uchime3_denovo ${outdir}/filtered/${sampleid}_denoised.fasta --sizein --sizeout --fasta_width 0 --nonchimeras ${outdir}/filtered/${sampleid}_denoised_nonchimeras.fasta

        ## starting from beginning, merge reads instead of concatenating them
        vsearch --fastq_mergepairs ${file} --reverse ${file/R1/R2} --fastq_allowmergestagger --threads ${nthreads} --fastaout - |

            ## fix headers to include sequence itself and sample ID
            awk -v sample=${sampleid} '$0 !~ ">" && length($0) > 0 {print ">"$0";sample="sample"\n"$0}' |

                ## derep sequences
                vsearch --derep_fulllength - --sizeout --fasta_width 0 --output ${outdir}/shortreads/merged/${sampleid}_derep.fasta
	fi

done

## derep across entire dataset
cat ${outdir}/filtered/*_denoised_nonchimeras.fasta |
    vsearch --derep_fulllength - --xsize --fasta_width 0 --output ${outdir}/rep_set.fasta

## create OTU table
cat ${outdir}/filtered/*_denoised_nonchimeras.fasta |
    vsearch --search_exact - --threads ${nthreads} -db ${outdir}/rep_set.fasta --sizein --sizeout --fasta_width 0 --otutabout ${outdir}/otu_table.txt

biom summarize-table -i ${outdir}/otu_table.txt -o ${outdir}/otu_table_summary.txt



## derep across entire dataset
cat ${outdir}/shortreads/merged/*.fasta |
    vsearch --derep_fulllength - --xsize --fasta_width 0 --output ${outdir}/shortreads/rep_set.fasta

## create OTU table
cat ${outdir}/shortreads/merged/*.fasta |
    vsearch --search_exact - --threads ${nthreads} -db ${outdir}/shortreads/rep_set.fasta --sizein --sizeout --fasta_width 0 --otutabout ${outdir}/shortreads/otu_table.txt

biom summarize-table -i ${outdir}/shortreads/otu_table.txt -o ${outdir}/shortreads/otu_table_summary.txt

