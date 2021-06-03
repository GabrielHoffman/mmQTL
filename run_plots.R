# Gabriel Hoffman
#
# October 2, 2020
# Intersect mmQTL with CAUSALdb

# Alz fine map
# https://www.medrxiv.org/content/10.1101/2020.01.22.20018424v1.full.pdf

# ml udunits proj gdal geos pandoc git; cd build2/mmQTL/; git pull
# source("run_plots.R")

# Change filt path of images
# sed -i 's/mmqtl/brema\/v2\/figures/g' *.html

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

synLogin()

source("./make_plots.R")
source("./plot_genes.R")
})
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

# df_snp = fread( paste0(folder, 'variants_position_from_MMQTL_meta-eQTL.tsv.gz') )
# colnames(df_snp) = c("ID", "Chr", "Position", "REF", "ALT")
# setkey(df_snp, 'ID')

df_snp = fread(synGet('syn25592274')$path)
colnames(df_snp) = c("ID", "Chr", "Position", "REF", "ALT")
setkey(df_snp, 'ID')

# df_eqtl = fread( paste0(folder, 'merge_MMQTL_eQTL_signal.tsv.gz'))
# df_eqtl[,Chr := c()]
# df_eqtl[,Z_score_fixed:=c()]
# setkey(df_eqtl, "Variant")
# gc()
# df_eqtl[,p.value := 2*pnorm(Z_score_random, 0, 1, lower.tail=FALSE)]

df_eqtl = fread(synGet('syn25592272')$path)
df_eqtl[,z_score_fixed:=NULL]
df_eqtl[,z_score_random:=NULL]
df_eqtl[,p_fixed:=NULL]
df_eqtl[,log10.p.value := -log10(p_random),]
df_eqtl[,p_random:=NULL]
df_eqtl[,Dis_to_TSS:=NULL]
df_eqtl[,Chr:=NULL]
setkey(df_eqtl, "Variant")
gc()


# pnorm(2, 0, 1, lower.tail=FALSE, log.p=TRUE) / log(10) + log10(2)
# log10(2*pnorm(2, 0, 1, lower.tail=FALSE)) 
# create log10 p directly
# df_eqtl[,log10.p.value := -1*pnorm(Z_score_random, 0, 1, lower.tail=FALSE, log.p=TRUE) / log(10) - log10(2)]
# df_eqtl[,Z_score_random:=c()]
# gc()

# table(df_eqtl$Variant %in% df_snp$ID)
df_eqtl = merge(df_eqtl, df_snp, by.x="Variant", by.y="ID")
setkey(df_eqtl, Position)
gc()

df_finemap = fread(synGet('syn25592273')$path)
# df_finemap = fread( paste0(folder, 'merge_MMQTL_eGenes_fine-mapping.tsv.gz') )
# df_finemap = fread( paste0(folder, 'merge_MMQTL_eGenes_fine-mapping_random_effect.tsv.gz') )
df_finemap[,Chr := c()]
# colnames(df_finemap)[4] = "PIP"
colnames(df_finemap)[colnames(df_finemap)=="PP"] = "PIP"

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
	file = "/sc/arion/projects/roussp01a/mads_atacseq/data/cellSpecificGenes/lake_habib_darmanis_2018.tsv"
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
df_merge = merge(df_finemap, df_causaldb, by="Variant", all=TRUE)

# merge with GWAS descriptions
df_merge = merge(df_merge, causalDB[['meta_info']], by.x="meta_id", by.y="ID")
df_merge[,PIP.prod := PIP*FINEMAP]
df_merge = df_merge[order(PIP.prod, decreasing=TRUE),]

# df_merge$Symbol = mapIds(org.Hs.eg.db,
#                      keys=df_merge$Gene,
#                      column="SYMBOL",
#                      keytype="ENSEMBL")

geneUnique = unique(df_merge$Gene)
geneUnique = geneUnique[!is.na(geneUnique)]

geneInfo = select(EnsDb.Hsapiens.v75, keys= geneUnique, keytype="GENEID", column=c('GENENAME'))
colnames(geneInfo) = c("Gene", "Symbol")	
df_merge = merge(df_merge, geneInfo, by="Gene", all.x=TRUE)



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

asdfasdfasdf

setwd("v2")

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
make_plot( ensGene, ord=1, showConditional= TRUE,window=2e5 )
make_plot( ensGene, ord=2, showConditional= TRUE,window=2e5 )

ensGene = 'ENSG00000140564'

make_plot( ensGene, ord=1, window=2e5 )



# source("../make_plots.R")
# make_plot( ensGene, ord=1, window=2e5 )


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

# saveRDS(df_show, file="df_show.RDS")

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



nrow(df_show)
length(table(df_show$Trait))
length(table(df_show$Variant))

# 

system("rsync -avzP figures/*.pdf sklar1:/sc/arion/projects/CommonMind/hoffman/MMQTL/figures")

system("rsync -avzP co_finemap_*.html sklar1:/sc/arion/projects/CommonMind/hoffman/MMQTL/figures")

# I integrated your finemappwing with GWAS finemapping from CAUSALdb (http://mulinlab.org/causaldb) to identify gene/trait pairs that share candidate causal variants.  Here is a list of the top associations with neuropsych traits: https://hoffmg01.u.hpc.mssm.edu/mmqtl/causal_coloc.html
# You can click the gene names for a detailed plot.  I am working on other traits too.  Are there any other GWAS finemapping datasets you want to add?

# Mechanism of variant
# TFBS
# splicing
# CpG

# Add eQTL and GWAS prob

####################
# FURIN pleiotropy #
####################

# df_show = readRDS("df_show.RDS")

df = df_show[Symbol=="FURIN",]
setorder(df, prob.coloc)
df = df[grep("predict|left|HURT", Trait,invert=TRUE),]
df$Trait = factor(df$Trait, df$Trait)

fig = ggplot(df, aes(Trait, prob.coloc)) + theme_bw() + geom_bar(stat="identity", fill="navy") + coord_flip() + scale_y_continuous(expand=c(0,0), limits=c(0,1)) + ylab("CLPP") + theme(aspect.ratio=1)
ggsave(file="FURIN_plieotropy.pdf", fig, height=5, width=5)








############################
# Summarize co-finemapping #
############################

# get MeSH_term
df2 = merge(df_show, 
	unique(causalDB[['meta_info']][,c('Trait', 'MeSH_term')]), 
	by='Trait')

tab = unique(df2[,data.frame(Count = length(unique(Gene))), by=c('CATEGORY', 'MeSH_term')])
tab$CATEGORY = factor(tab$CATEGORY)

setorder(tab, -CATEGORY, -Count)
tab$order = 1:nrow(tab)

# Show top result in each catagory
tab$show = TRUE
tab$show[-1] = sapply(2:nrow(tab), function(i) tab$CATEGORY[i] != tab$CATEGORY[i-1])

ylim = max(tab$Count)*1.05
fig = ggplot(tab, aes(order, Count, color=as.character(CATEGORY))) + geom_point() + theme_classic() + theme(aspect.ratio=2, plot.title = element_text(hjust = 0.5), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + ylab("# genes") + scale_x_continuous(expand=c(.01, 0)) + scale_y_continuous(expand=c(0,0), limits=c(0, ylim)) + geom_label_repel(data=tab[show==TRUE,], aes(order, Count, label=MeSH_term), force=4, nudge_x=10, nudge_y=10, size=2) + coord_flip()
ggsave(file = "figures/eQTL_GWAS_summary.pdf", fig, height=8, width=7)




# focus on traits
##################

keepTraits = c(
Neurotocism = "D000075384",
Alzheimers = 'D000544',
Schizophrenia = "D012559",
'Scz/Bipilar' = "D001523",
'Age first had sexual intercourse' = 'D003075',
'Worrier / anxious feelings' = 'D004644',
'Numeric memory test' = 'D000073219',
Intellegance = 'D007360',
Education = "D004522",
'Amyotrophic lateral sclerosis' = 'D000690',
'Multiple sclerosis' = 'D009103',
Alcohol = 'D000428',
'Anorexia nervosa' = 'D000856',
Smoking = 'D012907',
Bipolar = "D001714",
'Cognitive performance' = 'D003071',
'Depressive affect subcluster' = 'D000341',
Depression = 'D003863',
'Risk-Taking' = 'D012309',
Irritability = 'D007508',
'Mood swings' = 'D019964',
'Number of sexual partners' = 'D012725',
'Automobile speeding propensity' = 'D001334',
'Time spent watching television' = 'D013690',
'Coffee intake' = 'D003069',
Chronotype = 'D012890'
)



outfile = "brain_related_traits.html"

df_html[MeSH_ID %in% keepTraits,c('Category', 'Trait','Sample size', 'Author', 'PMID', 'Year', 'Gene','Symbol', 'Variant','GWAS prob', 'eQTL prob', 'Shared causal prob', 'Shared prob', 'Order', 'Cell type')] %>%
 	kable("html", escape = FALSE) %>%
  	kable_styling(full_width = FALSE, bootstrap_options = c("hover", "condensed")) %>%
  	save_kable(outfile)




tab = df2[MeSH_ID %in% keepTraits,]

tab$MeSH_term[tab$Trait=='Miserableness (MIS)'] = "Emotions"
tab$MeSH_term[tab$Trait=='Schizophrenia/Bipolar disorder'] = "Schizophrenia + Bipolar"
tab$MeSH_term[tab$Trait=='Schizophrenia'] = "Schizophrenia + Bipolar"
tab$MeSH_term[tab$Trait=='Bipolar Disorder'] = "Schizophrenia + Bipolar"
tab$MeSH_term[tab$Trait=='Bipolar disorder'] = "Schizophrenia + Bipolar"
tab$MeSH_term[tab$Trait=='Age first had sexual intercourse'] = "Age first had sexual intercourse"
tab$MeSH_term[grep('smoke',tab$Trait)] = "Smoking"
tab$MeSH_term[grep('elevision',tab$Trait)] = "Time spent watching television"
tab$MeSH_term[grep('speeding',tab$Trait)] = "Risk-Taking"
tab$MeSH_term[grep('Mood',tab$Trait)] = "Mood swings"
tab$MeSH_term[grep('Coffee',tab$Trait)] = "Coffee intake"

# tab[grep('isk-',Trait),] 
# tab[grep('Depression',MeSH_term),] 
tab[,Label := MeSH_term]
tab[,MeSH_term := c()]


tab[Label=='Affective Disorders, Psychotic',Label:='Depressive affect subcluster']



write.table( tab, file="Brain_eQTL_GWAS_focus.tsv", sep='\t', quote=FALSE)
write.table( df2, file="All_eQTL_GWAS.tsv", sep='\t', quote=FALSE)

tab_temp = read.delim("Brain_eQTL_GWAS_focus.tsv", sep="\t", header=TRUE)

# count brain related traits
nrow(tab_temp)
length(table(tab_temp$Label))
length(table(tab_temp$Variant))
length(table(tab_temp$Gene))

# count all traits
df2 =  read.delim("All_eQTL_GWAS.tsv.gz", sep="\t", header=TRUE)
nrow(df2)
length(table(df2$Trait))
length(table(df2$Variant))
length(table(df2$Gene))






tab2 = tab[,data.frame(Count = length(unique(Gene))), by='Label'] %>% unique

setorder(tab2, Count)
tab2$Label = factor(tab2$Label, tab2$Label)

ylim = max(tab2$Count)*1.05
fig = ggplot(tab2, aes(Label, Count, label=Count)) + geom_bar(stat="identity", fill="navy") + coord_flip() + theme_bw() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none") + scale_y_continuous(expand=c(0,0), limits=c(0, ylim)) + ylab("Genes with CLPP > 0.01") + xlab('') + geom_text(aes(y=Count + .6))

ggsave(file = "figures/Brain_eQTL_GWAS_focus_counts.pdf", fig, height=7, width=7)



tab3 = tab[MeSH_term %in% tab2$MeSH_term, data.frame(Summary = paste(sort(unique(Trait)), collapse=',')),by='MeSH_term']


# subset of candidate causal variants
#################
library(cowplot)
tab4 = tab[grep("Schiz|Multip|Alzh|Affect|Depr|Anor|Lateral", Label),]
setorder(tab4, Label, -prob.coloc)

tab4$Trait = factor(tab4$Trait, rev(sort(unique(tab4$Trait))))

figList = lapply( unique(tab4$Label), function(x){

	df = tab4[Label==x,]
	df[,Show := paste(Symbol, Variant, sep=' - ')] 
	df = df[,c('Show','prob.coloc', 'eQTL_order')][,.SD[which.max(prob.coloc)],by=Show]
	df$Show = factor(df$Show, rev(df$Show))

	fig = ggplot(df, aes(Show, prob.coloc, fill=as.character(eQTL_order))) + geom_bar(stat="identity") + theme_classic() + coord_flip() + scale_y_continuous(expand=c(0,0), limits=c(0, 1)) + scale_x_discrete(expand=c(0,0)) + ylab("CLPP") + ggtitle(x) + xlab('') + theme( legend.position="none", plot.title = element_text(hjust = 0.5))

	if( x != unique(tab4$Label)[length(unique(tab4$Label))] ){
			fig = fig + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
	}
	fig
})


fig = plot_grid( plotlist=figList, ncol=1, align="v", rel_heights=sapply(figList, function(x) nrow(x$data)))

ggsave(file = "figures/identify_Genes.pdf", fig, width=6, height=15)


df = tab4[Label=='Schizophrenia + Bipolar',]
df[,Show := paste(Symbol, Variant, sep=' - ')] 
df$Trait = gsub("Schizophrenia", "SCZ", df$Trait)
df$Trait = gsub("Bipolar disorder", "BD", df$Trait)
df2 = df[,c('Show','prob.coloc', 'eQTL_order')][,.SD[which.max(prob.coloc)],by=Show]
df2$Show = factor(df2$Show, rev(df2$Show))
df$Show = factor(df$Show, levels(df2$Show))

fig = ggplot(df, aes(Trait, Show, fill=prob.coloc)) + geom_tile() + theme_classic() + scale_fill_gradient(name="CLLP", low="lightgoldenrodyellow", high="red", limits=c(0, 1)) + ylab("") + theme(aspect.ratio=20/3)


ggsave(file = "figures/SCZ_BD_heatmap.pdf", fig, width=6, height=7)






tab2 = tab[,data.frame(Count = length(unique(Gene)), MeSH_term), by='MeSH_ID']
# tab2 = unique(tab[,data.frame(Count = length(unique(Gene)), Author=`Consortium/author`, N=Sample_size, Trait=Trait), by='MeSH_term'])

tab2$Trait = factor(tab2$Trait, unique(tab2$Trait[order(tab2$MeSH_term)]))

setorder(tab2, MeSH_ID)

ylim = max(tab2$Count)*1.05
ggplot(tab2, aes(Trait, Count, fill=MeSH_term)) + geom_bar(stat="identity") + coord_flip() + theme_classic() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="right") + scale_y_continuous(expand=c(0,0), limits=c(0, ylim)) + ylab("Genes with CLPP > 0.01")





# # keep
# idx = grep("Schiz|Alz|scler|eurot|ognativ|euromy|Alcohol|ducation|ntellig|sexual|memory|orrier|xiet", tab$Trait)
# tab2 = unique(tab[idx,])


# # remove
# idx2 = grep("doctor|Leisure", tab2$Trait, invert=TRUE)
# tab2 = tab2[idx2,]

# idx2 = grep("21926974|23453885|27046643", tab2$PMID, invert=TRUE)
# tab2 = tab2[idx2,]

# idx2
tab2 = tab2[idx2,] = grep("317694|26621|78308|264646|293723", tab2$N, invert=TRUE)


tab2 = tab[MeSH_ID %in% keepTraits,] 

tab2 = tab2[Count > 1,]

setorder(tab2, Trait)

ylim = max(tab2$Count)*1.05
ggplot(tab2, aes(paste(Trait, Author, N), Count, fill=MeSH_ID)) + geom_bar(stat="identity") + coord_flip() + theme_classic() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position="none") + scale_y_continuous(expand=c(0,0), limits=c(0, ylim)) 





idx = grep("Schiz", tab$Trait)

tab2 = tab[idx,]



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



# Load THOC7/FURIN example
####################

#  /sc/arion/projects/CommonMind/zengb02/finemapping_result_to_Gabriel.tar.gz

gene = "THOC7"

# folder = paste0("/Users/gabrielhoffman/Downloads/p_value_for_Gabriel/", gene,'/')
folder = "/Users/gabrielhoffman/Dropbox/projects/mmQTL/v2/from_Biao/Prepare_input_file_for_Gabriel_to_mmQTL_plot/Figure1/Figure1D/"

df_window = lapply( c('1', '7', '13', 'meta'), function(suffix){

	# For THOC
	file = paste0(folder,ifelse( suffix == "meta", "p_value_for_THOC7_meta", paste0("p_value_for_THOC7_", suffix,"_GTEx")))
	file2 = paste0(folder, ifelse(suffix == "meta", "finemapping_result_in_", "finemapping_result_in_GTEx_"), suffix)

	# if( suffix == "meta"){
	# 	file = paste0(folder,"eQTL_signal_in_meta")
	# 	file2 = paste0(folder,"finemapping_result_in_meta")
	# }else{
	# 	file = paste0(folder,"eQTL_signal_in_GTEX_", suffix)	
	# 	file2 = paste0(folder,"finemapping_result_in_GTEx_", suffix)		
	# }

	df_ex = fread( file )
	if( ncol(df_ex) == 4) df_ex = df_ex[,-1]
	df_ex$Description = suffix
	colnames(df_ex) = c("Variant", "beta", "p.value", "Description")
	df_ex = merge(df_ex, df_snp, by.x="Variant", by.y="ID")

	df_candSet = fread( file2, header=FALSE )
	colnames(df_candSet) = c('Variant')
	df_ex$inCandSet = "no"
	df_ex[Variant %in% df_candSet$Variant,inCandSet:="yes"]
	df_ex
})
df_window = do.call(rbind, df_window)

df_window[,min(p.value), by="Description"]

df_window[,sum(inCandSet=="yes"),by="Description"]



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

	ggplot(df_window[Description == dscr,], aes(Position, -log10(p.value), color=inCandSet)) + geom_point(size=1) + ggtitle("THOC7") + ylab(bquote(-log[10]~P)) + scale_x_continuous(expand=c(0,0), limits=c(start(wh), end(wh))) + scale_y_continuous(expand=c(0,0), limits=c(0, ymax*1.05), sec.axis = sec_axis(~./(ymax.rate/100), name = "Recombination\nrate [cM/Mb]")) + geom_line(aes(y=rate*(ymax.rate/100)), color="dodgerblue", size=.7) + theme_bw(8) + theme(aspect.ratio=.5, plot.title = element_text(hjust = 0.5)) +  scale_color_manual(values = c("black", "red")) #+ annotate("text", x=end(wh)- width(wh) / 5.5, y=ymax*0.9, label=paste("95% candidate set:", count), size=3)
})

library(EnsDb.Hsapiens.v75)
fig_genebody = plotEnsGenes_gg( EnsDb.Hsapiens.v75, start(wh), end(wh), seqnames(wh), splice_variants=FALSE) 

fig_main = make_plot( "ENSG00000163634", ord=1, window=window )

# fig_main = make_plot( "ENSG00000140564", ord=1, window=window )





# Combine plots
fig = tracks( "GTEx: 1" 	= figList[[1]],
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
 		title=gene) 
fig@mutable['Genes'] = FALSE

ggsave(file = "figures/example_THOC7.pdf", fig, width=6)
# ggsave(file = "figures/example_FURIN.pdf", fig, width=6)


# THOC7
# http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr3%3A63842598%2D63842660&hgsid=912394169_cFHqfAeeMj1FiG8bUqgEcpofXATB

APHB1 - rs117618017

# ZNF823 - rs72986630
# http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A11738636%2D11739205&hgsid=912424373_ZF4DWaMyOsuZrocMHs0reYmQysir

# disrupts signal in fetal brain
# https://hb.flatironinstitute.org/deepsea/jobs/32f5cf41-f21c-4f80-833c-6778a706318e

# https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr19%3A11731351%2D11746701&hgsid=912457981_vRgc0PjNmfRa16szbEvzWRRg3d0S

# Distrupt high information site in REST (aka NRSF) binding sites
# Regulome db

# CACHD1
# https://www.jneurosci.org/content/38/43/9186

# PPM1F-AS1
# https://doi.org/10.1016/j.abb.2018.01.001
# http://dx.doi.org/10.1016/j.biopsych.2017.08.013


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









