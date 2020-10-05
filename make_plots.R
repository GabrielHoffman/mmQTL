# Gabriel Hoffman
#
# October 2, 2020
# Intersect mmQTL with CAUSALdb

# Plot results
##############

make_plot = function( ensGene, window = 5e5, ord=1, non_coding=FALSE ){
	# FURIN
	# ZNF823
	# CD33
	# THOC7

	# ensGene = mapIds(org.Hs.eg.db,
	#          keys='ACE',
	#          column="ENSEMBL",
	#          keytype="SYMBOL")


	# ord = 1
	# window = 5e5

	geneSymbol = tryCatch(
			mapIds(org.Hs.eg.db,
	                     keys=ensGene,
	                     column="SYMBOL",
	                     keytype="ENSEMBL")
			, error = function(e) '')

	# get location of top fine-mapped eQTL SNP	
	df_fm = df_finemap[(Gene == ensGene) & (eQTL_order == ord),]
	pos_center = df_fm[which.max(PIP),Position]
	wh = GRanges(paste0("chr", df_fm$Chr[1]), IRanges(pos_center - window, pos_center + window))

	# eQTL
	df_window = df_eqtl[(Gene == ensGene) & (eQTL_order == ord),]

	# recombination rate
	df_recomb_sub = df_recomb[(position > start(wh)) & (position < end(wh)) & (chromosome == df_window$Chr[1]),]

	# Add recombination rate to first plot
	# linear interpolation evaluate at points in resComposite
	app = approx( df_recomb_sub$position, df_recomb_sub$recomb_rate, df_window$Position )
	df_window$rate = pmin(app$y, 100) # cap local rate at 100

	ymax = max(df_window$log10.p.value)
	ymax.rate = ymax
	# , color=as.character(eQTL_order)
	fig_eqtl = ggplot(df_window, aes(Position, log10.p.value)) + geom_point(size=.4) + ggtitle("FINEMAP - gene") + ylab(bquote(-log[10]~P)) + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0), limits=c(0, ymax*1.05), sec.axis = sec_axis(~./(ymax.rate/100), name = "Recombination\nrate [cM/Mb]")) + geom_line(aes(y=rate*(ymax.rate/100)), color="dodgerblue", size=.7) 

	# eQTL fine mapping
	# df_fm = df_finemap[(Gene == ensGene) & (eQTL_order == ord),]
	#, color=as.character(eQTL_order)
	fig_finemap_gene = ggplot(df_fm, aes(Position, PIP, size=PIP^2)) + geom_point() + ggtitle("Finemap - gene") + ylab("Posterior") + ylim(0,1) + scale_x_continuous(expand=c(0,0)) + scale_size_continuous(limits = c(0, 1), range(0.1, 6))

	# CAUSALdb
	# c('OT294', 'OT344', 'OT311', 'GA553', 'OT310')
	df_pip = get_credible_set( causalDB, unique(df_show$meta_id), wh)
	setorder(df_pip, Trait, -FINEMAP)

	df_coloc = merge(df_fm, df_pip, by.x="Position", by.y="BP")

	if( nrow(df_coloc) > 0){
		fig_causal = ggplot(df_coloc[FINEMAP*PIP > 0.01,], aes(Position, FINEMAP, color=Trait, size=FINEMAP^2)) + geom_point() + ylim(0,1) + scale_x_continuous(expand=c(0,0)) + ylab("Posterior") + geom_text_repel(data=subset(df_coloc, !duplicated(Trait)), aes(label=Trait), size=3, force=10) + scale_size_continuous(limits = c(0, 1), range(0.1, 6))
	}else{
		fig_causal = ggplot() + scale_x_continuous(expand=c(0,0), limits=c(start(wh), end(wh)))
	}

	# plot coloc
	# df_coloc = merge(df_fm, df_pip, by.x="Position", by.y="BP")

	if( nrow(df_coloc) > 0){
		ymax = max(df_coloc[,FINEMAP*PIP])
		fig_coloc = ggplot(df_coloc, aes(Position, FINEMAP*PIP, color=Trait, size=(FINEMAP*PIP)^2)) + geom_point() + scale_x_continuous(expand=c(0,0)) + ylab("Posterior") + scale_size_continuous(limits = c(0, 1), range(0.1, 6)) + ylim(0,1) #+ geom_text_repel(data=subset(df_coloc, !duplicated(Trait)), aes(label=Trait), size=4, force=10) 
	}else{	
		fig_coloc = ggplot() + scale_x_continuous(expand=c(0,0), limits=c(start(wh), end(wh)))
	}

	# Gene models
	# fig_genebody = tryCatch(
	# 		autoplot(db.ensdb.v75, wh, padding=0, fill="navy", color="navy", rect.height=.1) + scale_x_continuous(expand=c(0,0)) 
	# 		, error = function(e) ggplot())

	# browser()
 	library(EnsDb.Hsapiens.v75)
     
	fig_genebody = plotEnsGenes_gg( EnsDb.Hsapiens.v75, start(wh), end(wh), seqnames(wh), splice_variants=FALSE, non_coding=non_coding) 
	
	# tracks(fig_eqtl, fig_gene)

	# Combine plots
	fig_track = tracks( "eQTL"		= fig_eqtl, 
		"eQTL\nFine map"	= fig_finemap_gene,
		"GWAS\nFine map" 	= fig_causal,
		"Colocalize"		= fig_coloc,
		"Genes"				= fig_genebody,
		xlim = wh,
		padding = unit(-.8, "lines"),
		label.bg.fill="navy", label.text.color="white",
		heights=c(1,.5,.5,.5,.3),
		theme = theme_bw(8) + theme(legend.position="none"),
			title=paste0( geneSymbol, ' (',ensGene,')') )
	fig_track@mutable['Genes'] = FALSE
	fig_track
}


  # tracks(fig_eqtl, fig_finemap_gene, xlim=wh)









