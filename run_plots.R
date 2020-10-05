# Gabriel Hoffman
#
# October 2, 2020
# Intersect mmQTL with CAUSALdb

# Alz fine map
# https://www.medrxiv.org/content/10.1101/2020.01.22.20018424v1.full.pdf

# ml udunits proj gdal geos

library(data.table)
library(ggplot2)
library(ggrepel)
library(ggbio)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(dplyr)
library(knitr)
library(kableExtra)
library(decorate)
library(MeSH.db)

source("./make_plots.R")
source("./plot_genes.R")
 # make_plot( ensGene, ord=ord ) 


# plotEnsGenes_gg( EnsDb.Hsapiens.v75, 18565680, 19565680, "chr6", splice_variants=FALSE, non_coding=FALSE) 

# Read eQTL data
################

if( system("echo $HOSTNAME", intern=TRUE) == "glia.local"){
	folder = "/Users/gabrielhoffman/Downloads/MMQTL_eQTL_signal_and_finemapping_result/"
}else{
	folder = "/sc/arion/projects/CommonMind/hoffman/MMQTL/MMQTL_eQTL_signal_and_finemapping_result/"
}

df_snp = fread( paste0(folder, 'variants_position_from_MMQTL_meta-eQTL.tsv.gz') )
colnames(df_snp) = c("ID", "Chr", "Position", "REF", "ALT")
setkey(df_snp, 'ID')

df_eqtl = fread( paste0(folder, 'merge_MMQTL_eQTL_signal.tsv.gz'))
df_eqtl[,Chr := c()]
df_eqtl[,Z_score_fixed:=c()]
setkey(df_eqtl, "Variant")
gc()
# df_eqtl[,p.value := 2*pnorm(Z_score_random, 0, 1, lower.tail=FALSE)]

# pnorm(2, 0, 1, lower.tail=FALSE, log.p=TRUE) / log(10) + log10(2)
# log10(2*pnorm(2, 0, 1, lower.tail=FALSE)) 
# create log10 p directly
df_eqtl[,log10.p.value := -1*pnorm(Z_score_random, 0, 1, lower.tail=FALSE, log.p=TRUE) / log(10) - log10(2)]
df_eqtl[,Z_score_random:=c()]
gc()

# table(df_eqtl$Variant %in% df_snp$ID)
df_eqtl = merge(df_eqtl, df_snp, by.x="Variant", by.y="ID")
setkey(df_eqtl, Position)
gc()

df_finemap = fread( paste0(folder, 'merge_MMQTL_eGenes_fine-mapping.tsv.gz') )
df_finemap[,Chr := c()]
colnames(df_finemap)[4] = "PIP"

# table(df_finemap$Variant %in% df_snp$ID)
df_finemap = merge(df_finemap, df_snp, by.x="Variant", by.y="ID")
setkey(df_finemap, Position)

#  format(object.size(df_eqtl), "Mb")





# read CAUSALdb 
###############
# /Users/gabrielhoffman/workspace/scripts/compositeQTL/causaDB.R

if( system("echo $HOSTNAME", intern=TRUE) == "glia.local"){
	file = "/Users/gabrielhoffman/work/eQTL/resources/causaldb/credible_set/credible_set/causalDB.RDS"
}else{
	file = "/sc/arion/projects/psychencode/resources/causaldb/credible_set/causalDB.RDS"
}

causalDB = readRDS( file )

# find_study( causalDB, "Schizo")
find_study = function(causalDB, phenotypeText){
	causalDB[['meta_info']][grep(phenotypeText, causalDB[['meta_info']][['MeSH_term']]),]
}

# find_study( causalDB, "Alz")
# find_study( causalDB, "Park")

# find_study( causalDB, "Schi")


# get_credible_set( causalDB, 'OT294')

# get_credible_set( causalDB, c('OT294', 'AT095'))

# get_credible_set( causalDB, 'OT294')
get_credible_set = function(causalDB, studyID, gr){
	df = lapply( studyID, function(id){
		df = causalDB[["credible_set"]][meta_id==id,]
		df$Trait = causalDB[['meta_info']][ID==id,Trait] #MeSH_term
		df
		})
	df = do.call("rbind", df)	

	# if GenomeRanges object is specified, only return results in that range
	if( ! missing(gr) ){
		seqlevelsStyle(gr) = "NCBI"
		df_sub = df[CHR == as.vector(seqnames(gr)),]
		df_sub = df_sub[(BP > start(gr)) & (BP < end(gr)),]

		return( df_sub )
	}

	return( df )
}

# read gene models
##################

if( system("echo $HOSTNAME", intern=TRUE) == "glia.local"){
	gff_collapse = "/Users/gabrielhoffman/work/eQTL/resources/gene_models/ensembl_collapse/EnsDb.Hsapiens.v75.merged.gff3"
}else{
	gff_collapse = "/sc/arion/projects/psychencode/resources/gene_models/ensembl_collapse/EnsDb.Hsapiens.v75.merged.gff3"
}

db.ensdb.v75 = makeTxDbFromGFF( gff_collapse )


# Read Recombination rate map
##############################

if( system("echo $HOSTNAME", intern=TRUE) == "glia.local"){
	file = "/Users/gabrielhoffman/work/eQTL/resources/genetic_map/from_locuszoom/recombRate.tsv.gz"
}else{
	file = "/sc/arion/projects/psychencode/resources/genetic_map/from_locuszoom/recombRate.tsv.gz"
}

df_recomb = fread(file)

# intersect MMQTL with CAUSALdb
###############################

# Find which genes to plot
df_causaldb = causalDB[["credible_set"]]
df_causaldb[,Variant:=paste0('rs', rsID)]

# merge eQTL finemaping with causalDB
df_merge = merge(df_finemap, df_causaldb, by="Variant")

# merge with GWAS descriptions
df_merge = merge(df_merge, causalDB[['meta_info']], by.x="meta_id", by.y="ID")
df_merge[,PIP.prod := PIP*FINEMAP]
df_merge = df_merge[order(PIP.prod, decreasing=TRUE),]

df_merge$Symbol = mapIds(org.Hs.eg.db,
                     keys=df_merge$Gene,
                     column="SYMBOL",
                     keytype="ENSEMBL")



dim(df_merge[PIP.prod > 0.05,])
# sort(unique(df_merge[PIP.prod > 0.1,]$Trait))


res <- select(MeSH.db, 
			keys 	= unique(df_merge$MeSH_ID), 
			columns = c("CATEGORY", "MESHID" ,"MESHTERM"),
	 		keytype = 'MESHID')
res$CATEGORY = factor(res$CATEGORY)
res = res[order(res$CATEGORY),] 
res = res[!duplicated(res$MESHID),]

df_merge[(!df_merge$MeSH_ID %in% res$MESHID),]

df_merge = merge(df_merge, res, by.x="MeSH_ID", by.y='MESHID', all.x=TRUE)
df_merge$CATEGORY = as.character(df_merge$CATEGORY)
df_merge$CATEGORY[is.na(df_merge$CATEGORY)] = "Other"

# 177
# 201
# 215
# 273
# 438

# df_merge[273,]

# neuropsych traits
# phen = c('Alz',
# 'Cognitive',
# 'Drinks' ,
# 'Education',
# 'Intelligence',
# 'Lamb',
# 'sclerosis',
# 'euroticism',
# 'sexual',
# 'memory',
# 'Schiz',
# 'Worrier')
# reg = paste(phen, collapse="|")
# cutoff = 0.01
# outfile = "causal_coloc.html"

reg = ".*"
cutoff = 0.01
outfile = "causal_coloc_all.html"

df_show = df_merge[PIP.prod > cutoff,][grep(reg, Trait),c('CATEGORY', 'Trait','meta_id','Sample_size', 'Consortium/author', 'PMID', 'Year','PIP.prod', 'eQTL_order', 'Gene','Symbol')]

df_show = df_show[,.SD[which.max(PIP.prod),], by=c('Gene', 'Trait')]

df_show = df_show[order(CATEGORY, Trait, PIP.prod, decreasing=TRUE),][,c('CATEGORY', 'Trait','meta_id','Sample_size', 'Consortium/author', 'PMID', 'Year','PIP.prod', 'eQTL_order','Gene','Symbol')]

# Write images to pdf
for( ensGene in unique(df_show$Gene) ){
	message(ensGene)

	for(ord in df_show[Gene == ensGene,unique(eQTL_order)] ){

		# file = paste0("figures/", ensGene, "_", ord, ".pdf")
		file = paste0("/sc/arion/projects/CommonMind/hoffman/MMQTL/figures/", ensGene, "_", ord, ".pdf")		

		ggsave(file, make_plot( ensGene, ord=ord ) )
	}
} 

df_show[,url := paste0('https://hoffmg01.u.hpc.mssm.edu/mmqtl/',Gene, '_', eQTL_order,'.pdf')]
df_show[,url_pmid := paste0('https://pubmed.ncbi.nlm.nih.gov/', PMID)]
setorder(df_show, CATEGORY, Trait, -PIP.prod)

df_show %>% 
	mutate(Gene = cell_spec(Gene, "html", link = url, color="blue")) %>%
	mutate(url = c()) %>%
	mutate(PMID = cell_spec(PMID, "html", link = url_pmid, color="blue")) %>%
	mutate(url_pmid = c()) %>%
	mutate(meta_id = c()) %>%
	mutate('Coloc prob' = format(PIP.prod, digits=3)) %>%
	mutate(PIP.prod = c()) %>%
 	kable("html", escape = FALSE) %>%
  	kable_styling(full_width = FALSE, bootstrap_options = c("hover", "condensed")) %>%
  	save_kable(outfile)



system("rsync -avzP figures/*.pdf sklar1:/sc/arion/projects/CommonMind/hoffman/MMQTL/figures")

system("rsync -avzP causal_coloc*.html sklar1:/sc/arion/projects/CommonMind/hoffman/MMQTL/figures")

# https://hoffmg01.u.hpc.mssm.edu/mmqtl/causal_coloc.html
https://hoffmg01.u.hpc.mssm.edu/mmqtl/causal_coloc_all.html


# I integrated your finemappwing with GWAS finemapping from CAUSALdb (http://mulinlab.org/causaldb) to identify gene/trait pairs that share candidate causal variants.  Here is a list of the top associations with neuropsych traits: https://hoffmg01.u.hpc.mssm.edu/mmqtl/causal_coloc.html
# You can click the gene names for a detailed plot.  I am working on other traits too.  Are there any other GWAS finemapping datasets you want to add?


# Load THOC7 example
####################

folder = "/Users/gabrielhoffman/Downloads/p_value_for_Gabriel/"

df_ex1 = fread(paste0(folder, "p_value_1_GTEx"))
df_ex1$Description = "1"

df_ex7 = fread(paste0(folder, "p_value_7_GTEx"))
df_ex7$Description = "7"

df_ex13 = fread(paste0(folder, "p_value_13_GTEx"))
df_ex13$Description = "13"

df_ex_all = fread(paste0(folder, "p_value_Meta"))
df_ex_all$Description = "All"

df_window = rbind(df_ex1, df_ex7, df_ex13, df_ex_all )
colnames(df_window) = c("Chr", "Variant", "Position", "p.value", "Description")

window = 3e5
pos_center = df_window[which.min(p.value),Position]
wh = GRanges(paste0("chr", df_window$Chr[1]), IRanges(pos_center - window, pos_center + window))

# recombination rate
df_recomb_sub = df_recomb[(position > start(wh)) & (position < end(wh)) & (chromosome == df_window$Chr[1]),]

# Add recombination rate to first plot
# linear interpolation evaluate at points in resComposite
app = approx( df_recomb_sub$position, df_recomb_sub$recomb_rate, df_window$Position )
df_window$rate = pmin(app$y, 100) # cap local rate at 100

df_window[,Description:=factor(Description, rev(c("All", "13", "7", "1")))]

ymax = max(-log10(df_window$p.value))
ymax.rate = ymax


figList = lapply( levels(df_window$Description), function(dscr){

	ggplot(df_window[Description == dscr], aes(Position, -log10(p.value))) + geom_point(size=1) + ggtitle("THOC7") + ylab(bquote(-log[10]~P)) + scale_x_continuous(expand=c(0,0), limits=c(start(wh), end(wh))) + scale_y_continuous(expand=c(0,0), limits=c(0, ymax*1.05), sec.axis = sec_axis(~./(ymax.rate/100), name = "Recombination\nrate [cM/Mb]")) + geom_line(aes(y=rate*(ymax.rate/100)), color="dodgerblue", size=.7) + theme_bw(8) + theme(aspect.ratio=.5, plot.title = element_text(hjust = 0.5)) 
})

# fig_genebody = tryCatch(
# 		autoplot(db.ensdb.v75, wh, padding=0, fill="navy", color="navy", rect.height=.1) + scale_x_continuous(expand=c(0,0)) 
# 		, error = function(e) ggplot())
fig_genebody = plotEnsGenes_gg( EnsDb.Hsapiens.v75, start(wh), end(wh), seqnames(wh), splice_variants=FALSE) 



# Combine plots
fig_THOC7 = tracks( "GTEx: 1" 	= figList[[1]],
		"GTEx: 7" 	= figList[[2]],
		"GTEx: 13" 	= figList[[3]],
		"All" 		= figList[[4]],
		"Genes" 	= fig_genebody,
		padding = unit(-.8, "lines"),
		label.bg.fill="navy", label.text.color="white",
		heights=c(1,1,1,1,.4),
		theme = theme_bw(8) + theme(legend.position="none"),
 		title="THOC7") 
fig_THOC7@mutable['Genes'] = FALSE

ggsave(file = "figures/example_THOC7.pdf", fig_THOC7)


# THOC7
# http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr3%3A63842598%2D63842660&hgsid=912394169_cFHqfAeeMj1FiG8bUqgEcpofXATB

# ZNF823
# http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A11738636%2D11739205&hgsid=912424373_ZF4DWaMyOsuZrocMHs0reYmQysir




