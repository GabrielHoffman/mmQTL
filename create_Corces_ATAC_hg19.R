# Jaro Bendl and Gabriel Hoffman
# April 26, 2021
# Create GRanges list for snATAC-seq from Corces, et al, 2020, in hg19


# system("git pull")

library(GenomicRanges)
library(rtracklayer)

###
# CORCES :: Load all 24 clusters and assign them cluster (=cell type) names
PEAKFILES_CORCES1_ROOT = "/sc/arion/projects/roussp01a/jaro/data/creativity_gwas/input/ryan_corces/"
PEAKFILES_CORCES1 = list(
  "1_Isocortical_Excitatory"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster1.overlap.optimal.narrowPeak.gz"),
  "2_Striatal_Inhibitory_Major"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster2.overlap.optimal.narrowPeak.gz"),
  "3_Hippocampal_Excitatory"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster3.overlap.optimal.narrowPeak.gz"),
  "4_Hippocampal_Excitatory"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster4.overlap.optimal.narrowPeak.gz"),
  "5_Nigral_Neurons"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster5.overlap.optimal.narrowPeak.gz"),
  "6_Nigral_Neurons"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster6.overlap.optimal.narrowPeak.gz"),
  "7_Neurons_Unclassified"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster7.overlap.optimal.narrowPeak.gz"),
  "8_OPCs"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster8.overlap.optimal.narrowPeak.gz"),
  "9_OPCs"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster9.overlap.optimal.narrowPeak.gz"),
  "10_Nigral_OPCs"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster10.overlap.optimal.narrowPeak.gz"),
  "11_Isocortical_Inhibitory"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster11.overlap.optimal.narrowPeak.gz"),
  "12_Striatal_Inhibitory_Minor"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster12.overlap.optimal.narrowPeak.gz"),
  "13_Astrocytes_Unclassified"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster13.overlap.optimal.narrowPeak.gz"),
  "14_Nigral_Astrocytes"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster14.overlap.optimal.narrowPeak.gz"),
  "15_Isocortical_Astrocytes"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster15.overlap.optimal.narrowPeak.gz"),
  "16_Striatal_Astrocytes"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster16.overlap.optimal.narrowPeak.gz"),
  "17_Astrocytes_Unclassified"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster17.overlap.optimal.narrowPeak.gz"),
  "18_Potential_Doublets"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster18.overlap.optimal.narrowPeak.gz"),
  "19_Oligodendrocytes"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster19.overlap.optimal.narrowPeak.gz"),
  "20_Oligodendrocytes"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster20.overlap.optimal.narrowPeak.gz"),
  "21_Oligodendrocytes"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster21.overlap.optimal.narrowPeak.gz"),
  "22_Oligodendrocytes"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster22.overlap.optimal.narrowPeak.gz"),
  "23_Oligodendrocytes"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster23.overlap.optimal.narrowPeak.gz"),
  "24_Microglia"=file.path(PEAKFILES_CORCES1_ROOT, "Cluster24.overlap.optimal.narrowPeak.gz")
)
​
# chain for liftOver
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

PEAKS_GR = c()
for(file in names(PEAKFILES_CORCES1)) {
	fileDf = read.csv(PEAKFILES_CORCES1[[file]], sep="\t", skip=1, header=F, stringsAsFactors=F)[,1:4]
	rownames(fileDf) = paste0("Peak_corces1_", file , "_", 1:nrow(fileDf))
	fileDf$name = rownames(fileDf)
	colnames(fileDf) = c("seqnames", "start", "end", "name")
	if(!startsWith(fileDf$seqnames[1], prefix="chr")) { 
	fileDf$seqnames = paste0("chr", fileDf$seqnames)
	}

	fileDf = fileDf[fileDf$seqnames %in% c(paste0("chr", seq(1:22)), "chrX", "chrY", 1:22),]
	peakset = reduce(makeGRangesFromDataFrame(fileDf))

	# liftOver hg38 -> hg19
	seqlevelsStyle(gr) = "UCSC"  # necessary
	gr_hg19 = liftOver(gr, ch)
	class(gr_hg19)
	gr_hg19 = unlist(gr_hg19)

  	PEAKS_GR[[paste0("corcesOpt_", file)]] = gr_hg19
}
​
###
# CORCES :: Derived merged clusters for related cell types 
PEAKS_MERGED_GR = list()
PEAKS_MERGED_GR$Excitatory_neurons = GenomicRanges::union(PEAKS_GR$corcesOpt_1_Isocortical_Excitatory, GenomicRanges::union(PEAKS_GR$corcesOpt_3_Hippocampal_Excitatory, PEAKS_GR$corcesOpt_4_Hippocampal_Excitatory))
PEAKS_MERGED_GR$Inhibitory_neurons = GenomicRanges::union(PEAKS_GR$corcesOpt_2_Striatal_Inhibitory_Major, GenomicRanges::union(PEAKS_GR$corcesOpt_11_Isocortical_Inhibitory, PEAKS_GR$corcesOpt_12_Striatal_Inhibitory_Minor))
PEAKS_MERGED_GR$Nigral_neurons = GenomicRanges::union(PEAKS_GR$corcesOpt_5_Nigral_Neurons, PEAKS_GR$corcesOpt_6_Nigral_Neurons)
PEAKS_MERGED_GR$Neurons_unclassified = PEAKS_GR$corcesOpt_7_Neurons_Unclassified
PEAKS_MERGED_GR$OPCs = GenomicRanges::union(PEAKS_GR$corcesOpt_8_OPCs, GenomicRanges::union(PEAKS_GR$corcesOpt_9_OPCs, PEAKS_GR$corcesOpt_10_Nigral_OPCs))
PEAKS_MERGED_GR$Astrocytes = GenomicRanges::union(PEAKS_GR$corcesOpt_13_Astrocytes_Unclassified, GenomicRanges::union(PEAKS_GR$corcesOpt_14_Nigral_Astrocytes, GenomicRanges::union(PEAKS_GR$corcesOpt_15_Isocortical_Astrocytes, GenomicRanges::union(PEAKS_GR$corcesOpt_16_Striatal_Astrocytes, PEAKS_GR$corcesOpt_17_Astrocytes_Unclassified))))
PEAKS_MERGED_GR$Oligodendrocytes = GenomicRanges::union(PEAKS_GR$corcesOpt_19_Oligodendrocytes, GenomicRanges::union(PEAKS_GR$corcesOpt_20_Oligodendrocytes, GenomicRanges::union(PEAKS_GR$corcesOpt_21_Oligodendrocytes, GenomicRanges::union(PEAKS_GR$corcesOpt_22_Oligodendrocytes, PEAKS_GR$corcesOpt_23_Oligodendrocytes))))
PEAKS_MERGED_GR$Microglia = PEAKS_GR$corcesOpt_24_Microglia
​
saveRDS(PEAKS_MERGED_GR, file="Corces_hg19_merged.RDS")











