source ~/.bashrc

nthreads=4

indir=output/S1_20180718_15_28_51
outdir=output/S1_20180718_15_28_51

gg_reference_tree=/raid1/home/micro/mcmindsr/labhome/databases/gg_13_8_otus/trees/97_otus_unannotated.tree
gg_reference_seqs=/raid1/home/micro/mcmindsr/labhome/databases/gg_13_8_otus/rep_set/97_otus.fasta

source activate qiime1

parallel_assign_taxonomy_uclust.py -i ${outdir}/rep_set.fasta -o ${outdir}/otu_table_taxonomy -O ${nthreads}

mkdir ${outdir}/constrained_tree

vsearch --usearch_global ${outdir}/rep_set.fasta --db ${gg_reference_seqs} --id 0.80 --maxhits 2 --dbmatched ${outdir}/constrained_tree/gg_13_8_closedref97_matches.fasta

filter_tree.py -i ${gg_reference_tree} -f ${outdir}/constrained_tree/gg_13_8_closedref97_matches.fasta -o ${outdir}/constrained_tree/gg_13_8_97_otus_filtered.tre

perl /raid1/home/micro/mcmindsr/scripts/TreeToConstraints.pl < ${outdir}/constrained_tree/gg_13_8_97_otus_filtered.tre > ${outdir}/constrained_tree/gg_13_8_97_otus_filtered_constraints.txt

mafft --retree 1 --maxiterate 0 --nofft --parttree <(cat ${outdir}/rep_set.fasta ${outdir}/constrained_tree/gg_13_8_closedref97_matches.fasta) > ${outdir}/constrained_tree/aligned_ref_and_denoised_seqs.fasta

FastTree -constraints ${outdir}/constrained_tree/gg_13_8_97_otus_filtered_constraints.txt -nt < ${outdir}/constrained_tree/aligned_ref_and_denoised_seqs.fasta > ${outdir}/constrained_tree/gg_constrained_fastttree.tre

