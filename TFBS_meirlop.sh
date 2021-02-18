# Gabriel Hoffman
# Feb 18, 2021
#
# Analysis of TFBS in peaks that are positvely ore negatively 
# 	correlated with gene expression
# See software at https://github.com/npdeloss/meirlop

# Fits logistic regression model where the response is pressence/absense of motif in each sequence at p < x and covariates are score, and PC of kmer counts


cd /sc/arion/scratch/hoffmg01/meirlop

# donload BED file
# ml python/3.6.2
# synapse get syn24873889

# Run analysis

# ml anaconda3/2020.11
# alias conda=/hpc/packages/minerva-centos7/anaconda3/2020.11/bin/conda
# # conda create -n myenv bwa star
# conda activate myenv

# conda environment
ml anaconda3
source activate myenv


# Download motifs
JASPAR_MOTIFS_FILE_URL="http://jaspar2018.genereg.net/download/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
# User agent string, for spoofing the user agent to download the motifs file
USER_AGENT="Mozilla/5.0 (X11; Linux i686 (x86_64)) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/55.0.2883.75 Safari/537.36"

# Download JASPAR-formatted motifs file
wget -O motifs.txt -U "$USER_AGENT" "$JASPAR_MOTIFS_FILE_URL"

# FASTA=/sc/arion/projects/H_PBG/REFERENCES/GRCh38/FASTA/GRCh38.primary_assembly.genome.fa

# format the BED file
BED=ABC_PeakToGene_Spearman_Rho.bed
sed 's/[:blank:]+/\t/g' $BED | awk  -v OFS="\t" '{print $0, "+"}' | sed 's/ //g' > formated.bed

# R code to extract sequences from BED file and reference genome
# and create FASTA with scores for each sequence
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)

gr <- import("formated.bed")
seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,gr)
names(seq) = with(gr, paste0(seqnames, ':', start, '-', end, ' ', score))
writeXStringSet(seq,file="result.fa")

# Note that meirlop should take BED+FASTA file, but this crashes

# p < 1e-6 for motif match
meirlop \
	--jobs 48 \
	--fa result.fa \
	--html \
	--sortabs \
	--scan \
	--length \
	--pval 1e-6 \
	motifs.txt \
	meirlop_output_neg6/

# p < 1e-5 for motif match
meirlop \
	--jobs 48 \
	--fa result.fa \
	--html \
	--sortabs \
	--scan \
	--length \
	--pval 1e-5 \
	motifs.txt \
	meirlop_output_neg5/

# p < 1e-4 for motif match
meirlop \
	--jobs 48 \
	--fa result.fa \
	--html \
	--sortabs \
	--scan \
	--length \
	--pval 1e-4 \
	motifs.txt \
	meirlop_output_neg4/

# up load results to Synapse
synapse add --parentid syn24875664 --name meirlop_neg6_results.tsv meirlop_output_neg6/lr_results.tsv

synapse add --parentid syn24875664 --name meirlop_neg5_results.tsv meirlop_output_neg5/lr_results.tsv

synapse add --parentid syn24875664 --name meirlop_neg4_results.tsv meirlop_output_neg4/lr_results.tsv





########################
# R code to make plots #
########################

library(synapser)
library(data.table)
library(qvalue)
library(ggplot2)
library(cowplot)

synLogin()

file = synGet('syn24875668')$path # p < 1e-4
file = synGet('syn24875666')$path # p < 1e-5
file = synGet('syn24875670')$path # p < 1e-6

minMotifCount = 300

# Number of motifs to show in plot
nResults = 10

plot_enrichment = function(file, nResults, minMotifCount){

	df = fread( file )
	df$TFBS = gsub("^(.+) ", "", df$motif_id)

	df = df[num_peaks > minMotifCount,][order(coef),]

	df2 = data.frame(rbind(head(df, nResults), tail(df, nResults)))
	df2 = df2[order(df2$coef),]
	df2$TFBS = factor(df2$TFBS, df2$TFBS)

	lim = with(df2, max(abs(coef)))
	ggplot(df2, aes(TFBS,coef, fill=coef)) + geom_bar(stat="identity") + theme_classic(16) + coord_flip() + scale_fill_gradient2(name = "logOR", low="blue", mid="grey", high="red", limits=c(-lim, lim)) + ylab("log Odds Ratio") + xlab("TF motif") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1, legend.position="bottom")
}


file = synGet('syn24875668')$path # p < 1e-4
fig1 = plot_enrichment( file, nResults, minMotifCount) + ggtitle("p < 1e-4")

file = synGet('syn24875666')$path # p < 1e-5
fig2 = plot_enrichment( file, nResults, minMotifCount) + ggtitle("p < 1e-5")

file = synGet('syn24875670')$path # p < 1e-6
fig3 = plot_enrichment( file, nResults, minMotifCount) + ggtitle("p < 1e-6")

fig = plot_grid(fig1, fig2, fig3, nrow=1)

pdf("TFBS_enrichment.pdf", height=4, width=12)
fig
dev.off()


png("TFBS_enrichment.png", height=480, width=480*3)
fig
dev.off()








# ggplot(df2, aes(TFBS,coef/std_err, fill=log10(num_peaks))) + geom_bar(stat="identity") + theme_classic() + coord_flip() + scale_fill_gradient(low="grey", high="red")


# res = qvalue(df$pval)
 
# head(sort(res$qvalue))

# sum(res$qvalue < 0.2)

















