# Gabriel Hoffman
# Feb 2, 2021
#
# For eQTL and caQTL from the same set of SNPs
# 	evaluate how PIP scores from these assays relate



# relate eQTL and caQTL fine-mapping
# given PIP_e > x, what proprtion of PIP_ca > y

library(data.table)
library(ggplot2)
library(cowplot)
library(vcd)
library(synapser)

synLogin()

# caQTL fine-mapping results
df_caqtl = fread( synGet('syn24201357')$path, header=FALSE)
colnames(df_caqtl) = c("Chr", "Gene", "eQTL_order", "Variant", "PIP")
df_caqtl[,Chr := c()]
setkey(df_caqtl, 'Variant')
df_caqtl_best = df_caqtl[,.SD[which.max(PIP),],by="Variant"]

# eQTL fine-mapping results
df_eqtl = fread( synGet('syn24178479')$path )
colnames(df_eqtl)[5] = "PIP"
df_eqtl[,Chr := c()]
setkey(df_eqtl, 'Variant')
df_eqtl_best = df_eqtl[,.SD[which.max(PIP),],by="Variant"]

df_merge = merge(df_eqtl, df_caqtl_best, by="Variant", all=TRUE)
df_merge[is.na(PIP.x), PIP.x:=0]
df_merge[is.na(PIP.y), PIP.y:=0]
df_merge[,Gene.x:=c()]
df_merge[,Gene.y:=c()]
df_merge[,eQTL_order.x:=c()]
df_merge[,eQTL_order.y:=c()]
setkeyv(df_merge, c('PIP.x', 'PIP.y'))

s = seq(1e-3, .999, length.out=80)
grid = expand.grid(s,s)
nrow(grid)

df_count = lapply( seq_len(nrow(grid)), function(i){
	if( i %% round(nrow(grid)/30, 0) == 0) message(i)
	tab = df_merge[,table(PIP.x > grid$Var1[i], PIP.y > grid$Var2[i])]

	res = loddsratio(tab, log=TRUE)

	data.frame(	Var1 	= grid$Var1[i],
				Var2 	= grid$Var2[i], 
				a = tab[2,2] / rowSums(tab)[2],
				b = tab[2,2] / colSums(tab)[2],
				logOR 	= res$coefficients[1])
})
df = do.call(rbind, df_count)

fig1 = ggplot(df, aes(Var1, Var2, fill=logOR)) + geom_tile() + theme_classic() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position ="bottom") + scale_x_continuous(expand=c(0, 0), limits=c(0, 1)) + scale_y_continuous(expand=c(0, 0), limits=c(0, 1))+ scale_fill_gradient2( low="blue", mid="white", high="red") + xlab("eQTL PIP > x") + ylab("caQTL PIP > y") + ggtitle("log odds ratio")

lim = max(c(df$a, df$b))
fig2 = ggplot(df, aes(Var1, Var2, fill=a)) + geom_tile() + theme_classic() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position ="bottom") + scale_x_continuous(expand=c(0, 0), limits=c(0, 1)) + scale_y_continuous(expand=c(0, 0), limits=c(0, 1)) + scale_fill_gradient(name="Fraction of eQTL variants", low="white", high="red", limits=c(0, lim)) + xlab("eQTL PIP > x") + ylab("caQTL PIP > y") + ggtitle("Fraction of eQTL variants with co-incidence")

fig3 = ggplot(df, aes(Var1, Var2, fill=b)) + geom_tile() + theme_classic() + theme(aspect.ratio=1, plot.title = element_text(hjust = 0.5), legend.position ="bottom") + scale_x_continuous(expand=c(0, 0), limits=c(0, 1)) + scale_y_continuous(expand=c(0, 0), limits=c(0, 1))+ scale_fill_gradient(name="Fraction of caQTL variants", low="white", high="red", limits=c(0, lim)) + xlab("eQTL PIP > x") + ylab("caQTL PIP > y") + ggtitle("Fraction of caQTL variants with co-incidence")

pdf("QTL_coincidence.pdf", height=5, width=16)
plot_grid(fig1, fig2, fig3, nrow=1)
dev.off()


Test if these genes linked to peak with ABC are more likely to share a causal varint
































