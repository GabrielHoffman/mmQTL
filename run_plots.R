# Gabriel Hoffman
#
# October 2, 2020
# Intersect mmQTL with CAUSALdb

# Alz fine map
# https://www.medrxiv.org/content/10.1101/2020.01.22.20018424v1.full.pdf

# ml udunits proj gdal geos pandoc; cd build2/mmQTL/


library(data.table)
library(ggplot2)
library(ggrepel)
library(ggbio)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(dplyr)
library(stringr) 
library(knitr)
library(kableExtra)
library(decorate)
library(MeSH.db)

source("./make_plots.R")
source("./plot_genes.R")
 # make_plot( ensGene, ord=1 ) 


# plotEnsGenes_gg( EnsDb.Hsapiens.v75, 18565680, 19565680, "chr6", splice_variants=FALSE, non_coding=FALSE) 
# plotEnsGenes_gg( EnsDb.Hsapiens.v75, 193490008,194490008, "chr3", splice_variants=FALSE, non_coding=FALSE) 


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

# df_finemap = fread( paste0(folder, 'merge_MMQTL_eGenes_fine-mapping.tsv.gz') )
df_finemap = fread( paste0(folder, 'merge_MMQTL_eGenes_fine-mapping_random_effect.tsv.gz') )
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

# set keys for faster searching
setkey(causalDB[["credible_set"]], "meta_id")
setkey(causalDB[['meta_info']], "ID")

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
	# df = lapply( studyID, function(id){
	# 	df = causalDB[["credible_set"]][meta_id==id,]
	# 	df$Trait = causalDB[['meta_info']][ID==id,Trait] #MeSH_term
	# 	df
	# 	})
	# df = do.call("rbind", df)	

	# same as code above, but faster
	# df = causalDB[["credible_set"]][meta_id %in% studyID,]
	# df2 = causalDB[['meta_info']][ID %in% studyID,c('ID', 'Trait')]

	# use keys for speed
	df = causalDB[["credible_set"]][studyID,]
	df2 = causalDB[['meta_info']][studyID,c('ID', 'Trait')]
	df = merge(df, df2, by.x="meta_id", by.y="ID")

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


# single cell expression
########################

if( system("echo $HOSTNAME", intern=TRUE) == "glia.local"){
	file = "/Users/gabrielhoffman/work/eQTL/resources/singlecell/lake_habib_darmanis_2018.tsv"
}else{
	file = "/sc/hydra/projects/roussp01a/mads_atacseq/data/cellSpecificGenes/lake_habib_darmanis_2018.tsv"
}
df_sc = fread( file )

df_sc[df_sc$source=='darmanis',cluster:=clusterOriginal]
df_sc[df_sc$cluster=='',cluster:=clusterOriginal]
df_sc$source = str_to_title(df_sc$source)




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

# df_merge$Symbol = mapIds(org.Hs.eg.db,
#                      keys=df_merge$Gene,
#                      column="SYMBOL",
#                      keytype="ENSEMBL")

geneInfo = select(EnsDb.Hsapiens.v75, keys=df_merge$Gene, keytype="GENEID", column=c('GENENAME'))
colnames(geneInfo) = c("Gene", "Symbol")	
df_merge = merge(df_merge, geneInfo, by="Gene")



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
	# 'Amyotrophic',
	# Narcolepsy
	# Neuromyelitis
	# Alcohol
	# Anxiety
	# Bipolar
	# Depressi
	# smoke
	# CAT = "F"


# 'ducation',
# 'ntelligence',
# 'sclerosis',
# 'euroticism',
# 'sexual',
# 'memory',
# 'Schiz',
# 'Worrier')
# reg = paste(phen, collapse="|")
# cutoff = 0.01
# outfile = "causal_coloc.html"

# rs10015223


# Compute colocalization probability for each Gene, trait and eqtl order
df_show = merge(df_merge, 
			df_merge[, data.frame(prob.coloc = sum(PIP.prod)), by=c("Gene", 'meta_id', 'eQTL_order')], 
			by=c("Gene", 'meta_id', 'eQTL_order'))

reg = ".*"
cutoff = 0.01
# df_show[prob.coloc > cutoff,]
df_show = df_show[ PIP.prod> cutoff,][grep(reg, Trait),c('MeSH_ID','CATEGORY', 'Trait','meta_id','Sample_size', 'Consortium/author', 'PMID', 'Year', 'Variant','PIP','FINEMAP', 'PIP.prod', 'prob.coloc','eQTL_order', 'Gene','Symbol')]


df_show = df_show[,.SD[which.max(PIP.prod),], by=c('Gene', 'Trait')]

df_show = df_show[order(CATEGORY, Trait, PIP.prod, decreasing=TRUE),][,c('MeSH_ID','CATEGORY', 'Trait','meta_id','Sample_size', 'Consortium/author', 'PMID', 'Year','Variant','PIP','FINEMAP','PIP.prod', 'prob.coloc','eQTL_order','Gene','Symbol')]
df_show$Symbol[is.na(df_show$Symbol)] = ''

table(df_show$Symbol %in% df_sc$gene)

df_sc_unique = unique(df_sc[,c("gene", "source", "cluster")])
df_sc_unique[,cluster1 := paste0(cluster, ' (', source, ')')]
df_sc_unique = df_sc_unique[,data.table(clusters = paste(cluster1, collapse=', ')), by="gene"]

df_show = merge(df_show, df_sc_unique, by.x="Symbol", by.y="gene", all.x=TRUE)

df_show$clusters[is.na(df_show$clusters)] = ''

sdf

# Write images to pdf
for( ensGene in unique(df_show$Gene) ){
	message(ensGene)

	for(ord in df_show[Gene == ensGene,sort(unique(eQTL_order))] ){

		file = paste0(ensGene, "_", ord, ".pdf")
		file = paste0("figures/", file)
		# file = paste0("/sc/arion/projects/CommonMind/hoffman/MMQTL/figures/", file)		
		ggsave(file, make_plot( ensGene, ord=ord, window=2e5 ), width=6 )

		if( ord > 1 ){
			file = paste0(ensGene, "_", ord, "_showConditional.pdf")
			file = paste0("figures/", file)
			ggsave(file, make_plot( ensGene, ord=ord, window=2e5, showConditional=TRUE ), width=6) 
		}
	}
} 

# CACNA1C has two causal variants
ensGene = 'ENSG00000151067'
make_plot( ensGene, ord=1, showConditional= TRUE,window=2e6 )
make_plot( ensGene, ord=2, showConditional= TRUE,window=2e6 )



# rs240066
# ENSG00000224086 (an anti-sense RNA overlapping PPM1F) 
# ENSG00000100034 (PPM1F)
# ENSG00000100038 TOP3B

# df = df_finemap[Gene %in% c('ENSG00000224086', 'ENSG00000100034'),]
# setorder(df, -PIP)

# df_eqtl[Variant == 'rs240066',][order(log10.p.value),]


setorder(df_show, CATEGORY, Trait, -PIP.prod)

df_show[,url := paste0('https://hoffmg01.u.hpc.mssm.edu/mmqtl/',Gene, '_', eQTL_order,'.pdf')]
df_show[,url_pmid := paste0('https://pubmed.ncbi.nlm.nih.gov/', PMID)]
df_show[,url_rsid := paste0('https://www.ncbi.nlm.nih.gov/snp/', Variant)]
df_show[,url_mesh := paste0('http://id.nlm.nih.gov/mesh/', MeSH_ID)]
df_show[,url_symbol := paste0('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', Symbol)]
df_show[,url_gene := paste0('https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=', Gene)]

i = grep("Nealelab", df_show$PMID)
df_show$url_pmid[i] = "http://www.nealelab.is/uk-biobank"

# save(list=ls(), file="alldata.RDATA")

df_html = df_show %>% 
	mutate(Gene = cell_spec(Gene, "html", link = url_gene, color="blue", new_tab=TRUE)) %>%
	mutate(url_gene = c()) %>%
	mutate(PMID = cell_spec(PMID, "html", link = url_pmid, color="blue", new_tab=TRUE))%>%
	mutate(url_pmid = c()) %>%
	mutate(Variant = cell_spec(Variant, "html", link = url_rsid, color="blue", new_tab=TRUE)) %>%
	mutate(url_rsid = c()) %>%
	mutate(Symbol = cell_spec(Symbol, "html", link = url_symbol, color="blue", new_tab=TRUE)) %>%
	mutate(url_rsid = c()) %>%
	mutate(Trait = cell_spec(Trait, "html", link = url_mesh, color="blue", new_tab=TRUE)) %>%
	mutate(url_mesh = c()) %>%
	mutate('Sample size' = comma(as.integer(Sample_size), accuracy=1)) %>%
	rename(Author = `Consortium/author`) %>% 
	mutate(meta_id = c()) %>%
	rename(Category = CATEGORY) %>%
	mutate(`Shared causal prob` = cell_spec(format(PIP.prod, digits=3), "html", link = url, color="blue", new_tab=TRUE)) %>%
	mutate(PIP.prod = c())  %>%
	rename('Cell type' = clusters) %>%
	rename(Order = eQTL_order) %>%
	mutate('eQTL prob' = format(PIP, digits=3)) %>%
	mutate('GWAS prob' = format(FINEMAP, digits=3)) %>%
	mutate('Shared prob' = format(prob.coloc, digits=3)) 


for( CAT in unique(df_html$Category)){

	outfile = paste0("co_finemap_", CAT, ".html")

	df_html[Category==CAT,c('Category', 'Trait','Sample size', 'Author', 'PMID', 'Year', 'Gene','Symbol', 'Variant','GWAS prob', 'eQTL prob', 'Shared causal prob', 'Shared prob', 'Order', 'Cell type')] %>%
	 	kable("html", escape = FALSE) %>%
	  	kable_styling(full_width = FALSE, bootstrap_options = c("hover", "condensed")) %>%
	  	save_kable(outfile)
}




# saveRDS(df_show, file="df_show.RDS")


system("rsync -avzP figures/*.pdf sklar1:/sc/arion/projects/CommonMind/hoffman/MMQTL/figures")

system("rsync -avzP co_finemap_*.html sklar1:/sc/arion/projects/CommonMind/hoffman/MMQTL/figures")

# https://hoffmg01.u.hpc.mssm.edu/mmqtl/causal_coloc.html
https://hoffmg01.u.hpc.mssm.edu/mmqtl/causal_coloc_all.html


# I integrated your finemappwing with GWAS finemapping from CAUSALdb (http://mulinlab.org/causaldb) to identify gene/trait pairs that share candidate causal variants.  Here is a list of the top associations with neuropsych traits: https://hoffmg01.u.hpc.mssm.edu/mmqtl/causal_coloc.html
# You can click the gene names for a detailed plot.  I am working on other traits too.  Are there any other GWAS finemapping datasets you want to add?

Mechanism of variant
TFBS
splicing
CpG


Add eQTL and GWAS prob

# Load THOC7 example
####################

#  /sc/arion/projects/CommonMind/zengb02/finemapping_result_to_Gabriel.tar.gz

folder = "/Users/gabrielhoffman/Downloads/p_value_for_Gabriel/"

df_window = lapply( c('1', '7', '13', 'meta'), function(suffix){

	file = paste0(folder,ifelse( suffix == "meta", "p_value_Meta", paste0("p_value_", suffix,"_GTEx")))

	df_ex = fread( file )
	df_ex$Description = suffix
	colnames(df_ex) = c("Chr", "Variant", "Position", "p.value", "Description")

	file = paste0(folder, "finemapping_result_to_Gabriel/credible_set_GTEx_", suffix)
	df_candSet = fread( file )
	colnames(df_candSet) = c('Variant', "PIP1", "PIP2")
	df_ex$inCandSet = "no"
	df_ex[Variant %in% df_candSet$Variant,inCandSet:="yes"]
	df_ex
})
df_window = do.call(rbind, df_window)




window = 3e5
pos_center = df_window[which.min(p.value),Position]
wh = GRanges(paste0("chr", df_window$Chr[1]), IRanges(pos_center - window, pos_center + window))

# recombination rate
df_recomb_sub = df_recomb[(position > start(wh)) & (position < end(wh)) & (chromosome == df_window$Chr[1]),]

# Add recombination rate to first plot
# linear interpolation evaluate at points in resComposite
app = approx( df_recomb_sub$position, df_recomb_sub$recomb_rate, df_window$Position )
df_window$rate = pmin(app$y, 100) # cap local rate at 100

df_window[,Description:=factor(Description, rev(c("meta", "13", "7", "1")))]

ymax = max(-log10(df_window$p.value))
ymax.rate = ymax

figList = lapply( levels(df_window$Description), function(dscr){

	count = sum(df_window[Description == dscr,inCandSet=="yes"])

	ggplot(df_window[Description == dscr,], aes(Position, -log10(p.value), color=inCandSet)) + geom_point(size=1) + ggtitle("THOC7") + ylab(bquote(-log[10]~P)) + scale_x_continuous(expand=c(0,0), limits=c(start(wh), end(wh))) + scale_y_continuous(expand=c(0,0), limits=c(0, ymax*1.05), sec.axis = sec_axis(~./(ymax.rate/100), name = "Recombination\nrate [cM/Mb]")) + geom_line(aes(y=rate*(ymax.rate/100)), color="dodgerblue", size=.7) + theme_bw(8) + theme(aspect.ratio=.5, plot.title = element_text(hjust = 0.5)) +  scale_color_manual(values = c("black", "red")) + annotate("text", x=end(wh)- width(wh) / 5.5, y=ymax*0.9, label=paste("95% candidate set:", count), size=3)
})

library(EnsDb.Hsapiens.v75)
fig_genebody = plotEnsGenes_gg( EnsDb.Hsapiens.v75, start(wh), end(wh), seqnames(wh), splice_variants=FALSE) 

fig_main = make_plot( "ENSG00000163634", ord=1, window=window )


# Combine plots
fig_THOC7 = tracks( "GTEx: 1" 	= figList[[1]],
		"GTEx: 7" 	= figList[[2]],
		"GTEx: 13" 	= figList[[3]],
		"All" 		= figList[[4]],
		"eQTL\nFine map" 	= fig_main@plot[[2]],
		"GWAS\nFine map" 	= fig_main@plot[[3]],
		"Shared\ncandidates"= fig_main@plot[[4]],
		"Genes" 	= fig_genebody,
		xlim = wh,
		padding = unit(-.65, "lines"),
		label.bg.fill="navy", label.text.color="white",
		heights=c(1,1,1,1,1, 1,1,.65),
		theme = theme_bw(8) + theme(legend.position="none",panel.grid.minor = element_blank()),
 		title="THOC7") 
fig_THOC7@mutable['Genes'] = FALSE

ggsave(file = "figures/example_THOC7.pdf", fig_THOC7, width=6)


# THOC7
# http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr3%3A63842598%2D63842660&hgsid=912394169_cFHqfAeeMj1FiG8bUqgEcpofXATB

# ZNF823
# http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A11738636%2D11739205&hgsid=912424373_ZF4DWaMyOsuZrocMHs0reYmQysir

# disrupts signal in fetal brain
https://hb.flatironinstitute.org/deepsea/jobs/32f5cf41-f21c-4f80-833c-6778a706318e

https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A11731351%2D11746701&hgsid=912457981_vRgc0PjNmfRa16szbEvzWRRg3d0S

# Distrupt high information site in REST (aka NRSF) binding sites
# Regulome db

# CACHD1
https://www.jneurosci.org/content/38/43/9186

# PPM1F-AS1
https://doi.org/10.1016/j.abb.2018.01.001
http://dx.doi.org/10.1016/j.biopsych.2017.08.013


# ZNF823 brain ChIP-seq
library(rtracklayer)
library(ENCODExplorer)



query_results <- queryEncode(organism = "Homo sapiens", 
                      file_accession = "ENCFF625DED", file_format = "bigWig",
                      fixed = TRUE)

downloadEncode( query_results, dir="~/Downloads" )

library(rtracklayer)
library(ggplot2)
library(ggbio)
library(scales)

file = "/Users/gabrielhoffman/Downloads/ENCFF625DED.bigWig"

gr = import(file)

window = 1000
wh = GRanges("chr19", IRanges(11849736 - window, 11849736 + window))


gr2 = subsetByOverlaps(gr, wh, ignore.strand = TRUE)

df = as.data.frame(gr2)

ggplot(df, aes(start, score)) + geom_line() + scale_x_continuous(label = comma, expand=c(0,0), limits=c(start(wh), end(wh))) + theme_bw()









