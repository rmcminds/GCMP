source ~/.bashrc

source activate qiime1

# options: -d working directory, -t number of threads, -m mapping file, -k key to split table
nthreads=1
key=host_scientific_name
while getopts ":d:t:m:" opts
    do
    case $opts in
        d) workdir=$OPTARG ;;
        t) nthreads=$OPTARG ;;
        m) map=$OPTARG ;;
        k) key=$OPTARG ;;
    esac
done

outdir=output/${workdir}/

split_otu_table.py -i ${outdir}/shortreads/otu_table.txt -o ${outdir}/shortreads/per_host_otu_tables/ -m ${map} -f ${key}

mkdir ${outdir}/shortreads/per_host_otu_tables/mappings

mv ${outdir}/shortreads/per_host_otu_tables/*.txt ${outdir}/shortreads/per_host_otu_tables/mappings

for hostpath in ${outdir}/shortreads/per_host_otu_tables/*.biom ; do

    hostpath2=${hostpath/ /_}
    hostpath2=${hostpath2/per_host_otu_tables*${key}_/per_host_otu_tables\/}
    hostpath2=${hostpath2/_./.}
    hostpath2=${hostpath2/_./.}

    filter_otus_from_otu_table.py -i "${hostpath}" -o "${hostpath2%.*}_filt.biom" -n 1

    rm "${hostpath}"

    biom convert -i "${hostpath2%.*}_filt.biom" -o "${hostpath2%.*}.txt" --to-tsv --table-type 'OTU table' --header-key="taxonomy"

    rm "${hostpath2%.*}_filt.biom"

done
