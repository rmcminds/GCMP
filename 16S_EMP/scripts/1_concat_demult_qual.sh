source ~/.bashrc

# options: -d working directory, -t number of threads, -m mapping file, -i indir
nthreads=1
workdir=$(date +%C%y%m%d_%H_%M_%S)
while getopts ":d:i:m:t:" opts
    do
    case $opts in
        d) workdir=$OPTARG ;;
        i) indir=$OPTARG ;;
        m) sample2bc=$OPTARG ;;
        t) nthreads=$OPTARG ;;
    esac
done

outdir=output/${workdir}/

if [ ! -d "${outdir}" ]; then

    mkdir -p ${outdir}

    awk -F '\t' -v outdir=${outdir} 'NR>1 && $5=="GCMP1_Lane7_S7_L007" {print $1"\t"$2"\t""GCMP1" > outdir"/bcalls_GCMP1.txt"} NR>1 && $5=="GCMP2_Lane8_S8_L008" {print $1"\t"$2"\tGCMP2" > outdir"/bcalls_GCMP2.txt"}' ${sample2bc}

fi

if [ ! -d "${outdir}/demultiplexed/" ]; then

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

fi

if [ ! -d "${outdir}/filtered/" ]; then

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

                ## convert to fasta and construct headers as the sequence itself plus the sample ID
                awk -v sample=${sampleid} '(NR % 4) == 2 {print ">"$1$2";sample="sample"\n"$1$2}' |

                    ## derep sequences
                    vsearch --derep_fulllength - --sizeout --output - |

                        ## denoise sequences
                        vsearch --cluster_unoise - --sizein --sizeout --fasta_width 0 --minsize 2 --centroids - |

                            ## remove chimeras
                            vsearch --uchime3_denovo - --sizein --sizeout --fasta_width 0 --nonchimeras ${outdir}/filtered/${sampleid}_denoised_nonchimeras.fasta --chimeras ${outdir}/filtered/${sampleid}_denoised_chimeras.fasta

            ## starting from beginning, merge reads instead of concatenating them
            vsearch --fastq_mergepairs ${file} --reverse ${file/R1/R2} --fasta_width 0 --fastq_allowmergestagger --threads ${nthreads} --fastaout - |

                ## fix headers to include sequence itself and sample ID
                awk -v sample=${sampleid} '$0 !~ ">" && length($0) > 0 {print ">"$0";sample="sample"\n"$0}' |

                    ## derep sequences
                    vsearch --derep_fulllength - --sizeout --fasta_width 0 --output ${outdir}/shortreads/merged/${sampleid}_derep.fasta
        fi

    done

fi

## derep across entire dataset
cat ${outdir}/filtered/*_denoised_nonchimeras.fasta |
    vsearch --derep_fulllength - --xsize --fasta_width 0 --output ${outdir}/rep_set.fasta

## create OTU table
cat ${outdir}/filtered/*_denoised_nonchimeras.fasta |
    vsearch --search_exact - --threads ${nthreads} --db ${outdir}/rep_set.fasta --sizein --sizeout --fasta_width 0 --otutabout ${outdir}/otu_table.txt

biom summarize-table -i ${outdir}/otu_table.txt -o ${outdir}/otu_table_summary.txt



## derep across entire dataset
cat ${outdir}/shortreads/merged/*.fasta |
    vsearch --derep_fulllength - --xsize --fasta_width 0 --output ${outdir}/shortreads/rep_set.fasta

## create OTU table
cat ${outdir}/shortreads/merged/*.fasta |
    vsearch --search_exact - --threads ${nthreads} --minsize 2 -db ${outdir}/shortreads/rep_set.fasta --sizein --sizeout --fasta_width 0 --otutabout ${outdir}/shortreads/otu_table.txt

biom summarize-table -i ${outdir}/shortreads/otu_table.txt -o ${outdir}/shortreads/otu_table_summary.txt

