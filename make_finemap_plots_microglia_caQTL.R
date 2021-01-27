# Gabriel Hoffman
# January 15, 2021
#
# Dual fine mapping plots for microglia eQTL

suppressPackageStartupMessages({
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
library(synapser)
})

# ml proj gdal geos pandoc
# R
# source("/hpc/users/hoffmg01/build2/mmQTL/make_finemap_plots_microglia_caQTL.R")
  
src = '/hpc/users/hoffmg01/build2/mmQTL'

source(paste0(src, "/make_plots.R"))
source(paste0(src, "/plot_genes.R"))

synLogin()

# Read SNP locations
df_snp = fread( synGet('syn24201358')$path )
colnames(df_snp) = c("Chr", "Position0", "Position", 'ID')
df_snp[,Position0 := c()]
df_snp[,Chr:=gsub("^chr", '', Chr)]
setkey(df_snp, 'ID')
gc()

# Read xQTL p-values
df_eqtl = fread( synGet('syn24201360')$path, header=TRUE )
df_eqtl[,Chr:=gsub("^chr", '', Chr)]
colnames(df_eqtl)[colnames(df_eqtl) == 'caQTL_order'] = 'eQTL_order'
df_eqtl[,log10.p.value := -1*pnorm(Z_score_random, 0, 1, lower.tail=FALSE, log.p=TRUE) / log(10) - log10(2)]
df_eqtl[,Z_score_fixed:=c()]
df_eqtl[,Z_score_random:=c()]
colnames(df_eqtl)[colnames(df_eqtl)=="Peak"] = "Gene"
setkey(df_eqtl, "Variant")
gc()

# Include Position in df_eqtl
df_eqtl = merge(df_eqtl, df_snp[,c("Position", "ID")], by.x="Variant", by.y="ID")
setkey(df_eqtl, Position)
gc()

# Read fine-mapping results
df_finemap = fread( synGet('syn24201357')$path, header=FALSE)
colnames(df_finemap) = c("Chr", "Gene", "eQTL_order", "Variant", "PIP")
df_finemap[,Chr := c()]
setkey(df_finemap, 'Variant')

# Merg SNP locations and fine-mapping
df_finemap = merge(df_finemap, df_snp, by.x="Variant", by.y="ID")
setkey(df_finemap, Position)

ensGene = 'Peak_74506'

df_finemap[Gene == ensGene,]
df_eqtl[Gene == ensGene,]


table(unique(df_finemap$Gene) %in% unique(df_eqtl$Gene))




	




# read CAUSALdb 
###############
# /Users/gabrielhoffman/workspace/scripts/compositeQTL/causaDB.R

# This combines the results from CausalDB: http://mulinlab.tmu.edu.cn/causaldb/index.html
causalDB = readRDS( synGet('syn24178554')$path )

# set keys for faster searching
setkey(causalDB[["credible_set"]], "meta_id")
setkey(causalDB[['meta_info']], "ID")

find_study = function(causalDB, phenotypeText){
	causalDB[['meta_info']][grep(phenotypeText, causalDB[['meta_info']][['MeSH_term']]),]
}

get_credible_set = function(causalDB, studyID, gr){
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

# Look for any study with "Schizo" in the Trait 
# I want the study by Pardinas, et, al, 2018.
# This has the ID: OT294
res = find_study( causalDB, "Schizo")
head(res)

# Now get all the set of candidate causal variants for study OT294
df = get_credible_set( causalDB, 'OT294')

# There are a lot of columns, but this should be what you want.
# FINEMAP is the posterior inclusion probability
df[,c('meta_id','BP', 'rsID', 'FINEMAP', 'Trait')]


# read gene models
##################

db.ensdb.v75 = makeTxDbFromGFF( synGet('syn24194503')$path )

# Read Recombination rate map
##############################

df_recomb = fread( synGet('syn24194504')$path )

# single cell expression
########################

df_sc = fread( synGet('syn24194505')$path )

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

df_merge = df_merge[!is.na(PIP.prod),]

# df_merge$Symbol = mapIds(org.Hs.eg.db,
#                      keys=df_merge$Gene,
#                      column="SYMBOL",
#                      keytype="ENSEMBL")

# not needed for caQTLs
# geneInfo = select(EnsDb.Hsapiens.v75, keys=df_merge$Gene, keytype="GENEID", column=c('GENENAME'))
# colnames(geneInfo) = c("Gene", "Symbol")	
# df_merge = merge(df_merge, geneInfo, by="Gene")


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



# Compute colocalization probability for each Gene, trait and eqtl order
df_show = merge(df_merge, 
			df_merge[, data.frame(prob.coloc = sum(PIP.prod)), by=c("Gene", 'meta_id', 'eQTL_order')], 
			by=c("Gene", 'meta_id', 'eQTL_order'))

reg = ".*"
cutoff = 0.01
# df_show[prob.coloc > cutoff,]
df_show = df_show[ PIP.prod> cutoff,][grep(reg, Trait),c('MeSH_ID','CATEGORY', 'Trait','meta_id','Sample_size', 'Consortium/author', 'PMID', 'Year', 'Variant','PIP','FINEMAP', 'PIP.prod', 'prob.coloc','eQTL_order', 'Gene')]# ,'Symbol'


df_show = df_show[,.SD[which.max(PIP.prod),], by=c('Gene', 'Trait')]

df_show = df_show[order(CATEGORY, Trait, PIP.prod, decreasing=TRUE),][,c('MeSH_ID','CATEGORY', 'Trait','meta_id','Sample_size', 'Consortium/author', 'PMID', 'Year','Variant','PIP','FINEMAP','PIP.prod', 'prob.coloc','eQTL_order','Gene')]# ,'Symbol'
# df_show$Symbol[is.na(df_show$Symbol)] = ''

# table(df_show$Symbol %in% df_sc$gene)

# df_sc_unique = unique(df_sc[,c("gene", "source", "cluster")])
# df_sc_unique[,cluster1 := paste0(cluster, ' (', source, ')')]
# df_sc_unique = df_sc_unique[,data.table(clusters = paste(cluster1, collapse=', ')), by="gene"]

# df_show = merge(df_show, df_sc_unique, by.x="Symbol", by.y="gene", all.x=TRUE)

# df_show$clusters[is.na(df_show$clusters)] = ''


folder = "microglia/caqtl/figures/"
dir.create( folder, recursive=TRUE )

# Write images to pdf
for( ensGene in unique(df_show$Gene) ){
	message(ensGene)

	for(ord in df_show[Gene == ensGene,sort(unique(eQTL_order))] ){

		file = paste0(ensGene, "_", ord, ".pdf")
		file = paste0(folder, file)
		# file = paste0("/sc/arion/projects/CommonMind/hoffman/MMQTL/figures/", file)		
		ggsave(file, make_plot( ensGene, ord=ord, window=8e4 ), width=6 )

		if( ord > 1 ){
			file = paste0(ensGene, "_", ord, "_showConditional.pdf")
			file = paste0("figures/", file)
			ggsave(file, make_plot( ensGene, ord=ord, window=2e5, showConditional=TRUE ), width=6) 
		}
	}
} 



# ensGene = "ENSG00000128606"
# make_plot( ensGene, ord=1,window=1e6 )


# df_show[Gene == 'ENSG00000128606',]


# # CACNA1C has two causal variants
# ensGene = 'ENSG00000136717'
# make_plot( ensGene, ord=1, showConditional= TRUE,window=2e6 )
# make_plot( ensGene, ord=2, showConditional= TRUE,window=2e6 )


# Make HTML files
#################


setorder(df_show, CATEGORY, Trait, -PIP.prod)

df_show[,url := paste0('../figures/',Gene, '_', eQTL_order,'.pdf')]
df_show[,url_pmid := paste0('https://pubmed.ncbi.nlm.nih.gov/', PMID)]
df_show[,url_rsid := paste0('https://www.ncbi.nlm.nih.gov/snp/', Variant)]
df_show[,url_mesh := paste0('http://id.nlm.nih.gov/mesh/', MeSH_ID)]
# df_show[,url_symbol := paste0('https://www.genecards.org/cgi-bin/carddisp.pl?gene=', Symbol)]
df_show[,url_gene := paste0('https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=', Gene)]

i = grep("Nealelab", df_show$PMID)
df_show$url_pmid[i] = "http://www.nealelab.is/uk-biobank"

saveRDS(df_show, file="df_show.RDS")

# save(list=ls(), file="alldata.RDATA")

df_html = df_show %>% 
	mutate(Gene = cell_spec(Gene, "html", link = url_gene, color="blue", new_tab=TRUE)) %>%
	mutate(url_gene = c()) %>%
	mutate(PMID = cell_spec(PMID, "html", link = url_pmid, color="blue", new_tab=TRUE))%>%
	mutate(url_pmid = c()) %>%
	mutate(Variant = cell_spec(Variant, "html", link = url_rsid, color="blue", new_tab=TRUE)) %>%
	mutate(url_rsid = c()) %>%
	# mutate(Symbol = cell_spec(Symbol, "html", link = url_symbol, color="blue", new_tab=TRUE)) %>%
	mutate(url_rsid = c()) %>%
	mutate(Trait = cell_spec(Trait, "html", link = url_mesh, color="blue", new_tab=TRUE)) %>%
	mutate(url_mesh = c()) %>%
	mutate('Sample size' = comma(as.integer(Sample_size), accuracy=1)) %>%
	rename(Author = `Consortium/author`) %>% 
	mutate(meta_id = c()) %>%
	rename(Category = CATEGORY) %>%
	mutate(`Shared causal prob` = cell_spec(format(PIP.prod, digits=3), "html", link = url, color="blue", new_tab=TRUE)) %>%
	mutate(PIP.prod = c())  %>%
	# rename('Cell type' = clusters) %>%
	rename(Order = eQTL_order) %>%
	mutate('eQTL prob' = format(PIP, digits=3)) %>%
	mutate('GWAS prob' = format(FINEMAP, digits=3)) %>%
	mutate('Shared prob' = format(prob.coloc, digits=3)) 


dir.create( paste0(folder, "../html/"), recursive=TRUE )

for( CAT in unique(df_html$Category)){

	outfile = paste0(folder, "../html/co_finemap_", CAT, ".html")

	df_html[Category==CAT,c('Category', 'Trait','Sample size', 'Author', 'PMID', 'Year', 'Gene', 'Variant','GWAS prob', 'eQTL prob', 'Shared causal prob', 'Shared prob', 'Order')] %>%
	 	kable("html", escape = FALSE) %>%
	  	kable_styling(full_width = FALSE, bootstrap_options = c("hover", "condensed")) %>%
	  	save_kable(outfile)
}




























