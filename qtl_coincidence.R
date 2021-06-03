# Gabriel Hoffman
# Feb 2, 2021
#
# For eQTL and caQTL from the same set of SNPs
# 	evaluate how PIP scores from these assays relate



# relate eQTL and caQTL fine-mapping
# given PIP_e > x, what proprtion of PIP_ca > y

suppressPackageStartupMessages({
library(data.table)
library(ggplot2)
library(cowplot)
library(vcd)
library(GenomicRanges)
library(AnnotationHub)
library(ensembldb)
library(synapser)
})

synLogin()



# For ENSEMBL id ENSG00000279457.4, return ENSG00000279457
trim_ensembl_ids = function(x){
  gsub("(.*)\\.(.*)", "\\1", x) 
}

# Add column Symbol using column gene_id storing ENSEMBL id
getGeneSymbol = function( df, column="row.names"){

  if( ! is(df, "data.frame") ){
    df = as.data.frame(df)
  }

  # Load ENSEMBL v96 database
  # see https://www.bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/ensembldb.html#10_getting_or_building_ensdb_databasespackages
  ah <- AnnotationHub()
  ensdb = ah[["AH69187"]] # ENSEMBL v96

  if( column == "row.names"){    
    geneID = rownames(df)
  }else{
    geneID = df[[column]]
  }

  ensIDs = trim_ensembl_ids( geneID )

  geneSymbol = mapIds(ensdb, keys=ensIDs, keytype="GENEID", column="GENENAME")
  df$Symbol = geneSymbol

  entrez = mapIds(ensdb, keys=ensIDs, keytype="GENEID", column="ENTREZID")
  df$Entrez = entrez

  df_info = ensembldb::select(ensdb, keys=ensIDs, keytype="GENEID", column=c("SEQNAME", "GENESEQSTART", "GENESEQEND"))

  # df2 = merge( data.frame(GENEID = ensIDs, stringsAsFactors=FALSE), chroms, by='GENEID')

  idx = match(ensIDs, df_info$GENEID)
  # df$Chrom = c()
  df$Chrom = df_info$SEQNAME[idx]
  # df$Start = c()
  df$Start = df_info$GENESEQSTART[idx]
  # df$End = c()
  df$End = df_info$GENESEQEND[idx]

  df
}


# caQTL fine-mapping results
df_caqtl = fread( synGet('syn24911652')$path, header=FALSE) # syn24201357
colnames(df_caqtl) = c("Chr", "Peak", "eQTL_order", "Variant", "PIP")
df_caqtl[,Chr := c()]
setkey(df_caqtl, 'Variant')
df_caqtl_best = df_caqtl[,.SD[which.max(PIP),],by="Variant"]

# peak location
df_peak = fread( synGet('syn24201356')$path )
colnames(df_peak) = c("Chr", "start", "end", 'Peak')
gr_peak = with(df_peak, GRanges( Chr, IRanges(start, end, name=Peak)))


# eQTL fine-mapping results
df_eqtl = fread( synGet('syn24178479')$path )
colnames(df_eqtl)[5] = "PIP"
df_eqtl[,Chr := c()]
setkey(df_eqtl, 'Variant')
df_eqtl_best = df_eqtl[,.SD[which.max(PIP),],by="Variant"]

# location of genes
df_gene = getGeneSymbol( df_eqtl[,data.frame(Gene = unique(Gene))], "Gene")
gr_gene = with(df_gene, GRanges(Chrom, IRanges(Start, End, name=Gene)))

# ABC links
df_abc = readRDS( synGet('syn24346713')$path )$ALL_CPM1
# df_abc = df_abc[df_abc$ABC.Score > .1,]

# Merge eQTL and caQTL
df_merge = merge(df_eqtl, df_caqtl_best, by="Variant", all=TRUE)
df_merge[is.na(PIP.x), PIP.x:=0]
df_merge[is.na(PIP.y), PIP.y:=0]
df_merge[,eQTL_order.x:=c()]
df_merge[,eQTL_order.y:=c()]
df_merge[,combo := paste(Peak, Gene, sep='_') ]

# Indicate if gene-peak pair is an ABC link
df_merge[,isABC:=FALSE]
df_merge[combo %in% df_abc$combo,isABC:=TRUE]
setkeyv(df_merge, c('PIP.x', 'PIP.y'))


# Compute distance between genes and Peaks
idx = df_merge[,which(!is.na(Gene) & !is.na(Peak))]
dst = mclapply(idx, function(i){
	if( i %% 1000 == 0){
		frac = min(which(idx == i)) / length(idx)
		message(format(frac*100, digits=3), '%')
	}
	distance(gr_gene[df_merge$Gene[i]], gr_peak[df_merge$Peak[i]])
	}, mc.cores=6)
dst = unlist(dst)
df_merge[,distance :=0]
# df_merge[,distance :=NA]
df_merge$distance[idx] = dst


# Define functions so analysis can be repeated
compute_enrichment = function(df_merge, len){
	s = seq(1e-3, .999, length.out=len)
	grid = expand.grid(s,s)
	nrow(grid)

	df_count = mclapply( seq_len(nrow(grid)), function(i){
		if( i %% round(nrow(grid)/30, 0) == 0) message(i)
		tab = df_merge[,table(PIP.x > grid$Var1[i], PIP.y > grid$Var2[i])]

		if( length(tab) == 4){
			res = loddsratio(tab, log=TRUE)

			frac_eQTL 	= tab[2,2] / rowSums(tab)[2]
			frac_caQTL 	= tab[2,2] / colSums(tab)[2]
			logOR 		= res$coefficients[1]
			se 			= as.numeric(sqrt(vcov(res)))
		}else{
			frac_eQTL = NA
			frac_caQTL = NA
			logOR = NA
			se = NA
		}

		data.frame(	Var1 		= grid$Var1[i],
					Var2 		= grid$Var2[i], 
					frac_eQTL 	= frac_eQTL,
					frac_caQTL 	= frac_caQTL,
					logOR 		= logOR,
					se 			= se)
	}, mc.cores=6)
	data.table(do.call(rbind, df_count))
}

make_plots = function(df){
	lim = min(min(df$logOR), 0)
	fig1 = ggplot(df, aes(Var1, Var2, fill=logOR)) + geom_raster() + theme_classic() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position ="bottom") + scale_x_continuous(expand=c(0, 0), limits=c(0, 1)) + scale_y_continuous(expand=c(0, 0), limits=c(0, 1))+ scale_fill_gradient2( low="blue", mid="white", high="red", limits=c(lim, NA)) + xlab("eQTL PIP > x") + ylab("caQTL PIP > y") + ggtitle("log odds ratio")

	lim = max(c(df$frac_eQTL, df$frac_caQTL))
	fig2 = ggplot(df, aes(Var1, Var2, fill=frac_eQTL)) + geom_raster() + theme_classic() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position ="bottom") + scale_x_continuous(expand=c(0, 0), limits=c(0, 1)) + scale_y_continuous(expand=c(0, 0), limits=c(0, 1)) + scale_fill_gradient(name="Fraction of eQTL variants", low="white", high="red", limits=c(0, lim)) + xlab("eQTL PIP > x") + ylab("caQTL PIP > y") + ggtitle("Fraction of eQTL variants with co-incidence")

	fig3 = ggplot(df, aes(Var1, Var2, fill=frac_caQTL)) + geom_raster() + theme_classic() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position ="bottom") + scale_x_continuous(expand=c(0, 0), limits=c(0, 1)) + scale_y_continuous(expand=c(0, 0), limits=c(0, 1))+ scale_fill_gradient(name="Fraction of caQTL variants", low="white", high="red", limits=c(0, lim)) + xlab("eQTL PIP > x") + ylab("caQTL PIP > y") + ggtitle("Fraction of caQTL variants with co-incidence")
	plot_grid(fig1, fig2, fig3, nrow=1) 
}

max(df_merge$distance, na.rm=TRUE)


distCutoff = 1e4 
res_abc = compute_enrichment( df_merge[(isABC==TRUE) & (distance < distCutoff),], 15)
res_notabc = compute_enrichment( df_merge[(isABC==FALSE) & (distance < distCutoff),], 15)

# 


res_abc = compute_enrichment( df_merge[(isABC==TRUE),], 17)
res_notabc = compute_enrichment( df_merge[(isABC==FALSE),], 17)
# res_notabc = compute_enrichment( df_merge[(isABC==FALSE) & !is.na(distance),], 15)


fig1 = make_plots(res_notabc)
fig2 = make_plots(res_abc)

# Compare ABC versus others
res_combined = res_abc[,1:2]
res_combined$frac_eQTL = log(res_abc$frac_eQTL) - log(res_notabc$frac_eQTL) 
res_combined$frac_caQTL = log(res_abc$frac_caQTL) - log(res_notabc$frac_caQTL) 
res_combined$logOR = res_abc$logOR - res_notabc$logOR

fig3 = make_plots( res_combined )


pdf("QTL_coincidence.pdf", height=16, width=16)
plot_grid(fig1, fig2, fig3, ncol=1, labels = c("Non-ABC", "ABC", "Compare"))
dev.off()



# GEH May 5, 2021
#################


# Define functions so analysis can be repeated
compute_enrichment = function(df_merge, len){
	s = seq(1e-3, .999, length.out=len)

	df_count = lapply(s, function(value){
		# if( i %% round(nrow(grid)/30, 0) == 0) message(i)
		tab = df_merge[,table(PIP.x > value, PIP.y > value)]

		if( length(tab) == 4){
			res = loddsratio(tab, log=TRUE)

			frac_eQTL 	= tab[2,2] / rowSums(tab)[2]
			frac_caQTL 	= tab[2,2] / colSums(tab)[2]
			logOR 		= res$coefficients[1]
			se 			= as.numeric(sqrt(vcov(res)))
		}else{
			frac_eQTL = NA
			frac_caQTL = NA
			logOR = NA
			se = NA
		}

		data.frame(	Var1 		= value,
					Var2 		= value, 
					frac_eQTL 	= frac_eQTL,
					frac_caQTL 	= frac_caQTL,
					logOR 		= logOR,
					se 			= se)
	})
	data.table(do.call(rbind, df_count))
}

n_pts = 20
res_abc = compute_enrichment( df_merge[(isABC==TRUE),], n_pts)
res_notabc = compute_enrichment( df_merge[(isABC==FALSE),], n_pts)

df = lapply( unique(res_notabc$Var1), function(x){
	a = data.frame(isABC ='1', res_abc[Var1==x & Var2 ==x,])
	b = data.frame(isABC ='0', res_notabc[Var1==x & Var2 ==x,])
	rbind(a,b)
})
df = do.call(rbind, df)


ggplot(df, aes(Var1, logOR, fill=isABC)) + geom_ribbon(aes(ymin = logOR - 1.96*se, ymax = logOR + 1.96*se), width=.02, alpha=.2)  + geom_line(aes(, color=isABC)) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + theme_classic() + xlab("PIP cutoff") + scale_fill_manual(name="ABC link", values=c("grey30", 'red2')) + scale_color_manual(name="ABC link", values=c("grey30", 'red2')) + geom_hline(yintercept=0, color = "black", linetype="dashed") + ggtitle("Enrichment in overlap from eQTL and caQTL") + xlim(0, 1) + ylim(-4, 10)




# ggplot(df3, aes(Var1, logOR, color=isABC)) + geom_point() + geom_errorbar(aes(ymin = logOR - 1.96*se, ymax = logOR + 1.96*se), width=.02) + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5)) + theme_classic() + xlab("PIP cutoff") + scale_color_manual(name="ABC link",  values=c("grey30", 'red2')) + geom_hline(yintercept=0, color = "black", linetype="dashed")




























