
library(tidyverse)
library(stringr)
library(gtools)

# sum the number of Barcode, observing the same barcode multiple times doesn't contribute to the count

folder_mpra = "/sc/arion/projects/CommonMind/leed62/fine_mapping_SMI_MPRA/MPRA_results/"
# read in data
file = paste0(folder_mpra, "fine_mapping_final_labeled_counts.tsv")
df = read_tsv(file)
colnames(df)[4:9] = c('DNA1', 'DNA2', "DNA3", "RNA1", "RNA2", "RNA3")

# For each label, sum the number of barcodes observed 
# 	(i.e. with nonzero counts)
# Observing the same barcode multiple times doesn't
#	contribute to the count
df2 = df %>%
	filter(label != "no_BC")  %>%
	group_by(label) %>%
	summarize(DNA1 = sum(!is.na(DNA1)),
			DNA2 = sum(!is.na(DNA2)),
			DNA3 = sum(!is.na(DNA3)),
			RNA1 = sum(!is.na(RNA1)),
			RNA2 = sum(!is.na(RNA2)),
			RNA3 = sum(!is.na(RNA3)))

# get two control categories
df_controls_c1 = df2 %>% filter(grepl("^c1", label))
df_controls_c2 = df2 %>% filter(grepl("^c2", label))

# get negative sequences
df_sub_n = df2 %>% filter(grepl("^n", label))

# extract info from label
df_label_n = str_split(df_sub_n$label, '_|\\:\\:', simplify=TRUE)
colnames(df_label_n) = c("Category", "rsid", "ref", "alt", "variant", "position", 'genomeCoords')
df_label_n = data.frame(df_label_n)
df_label_n$label = df_sub_n$label

# get positive sequences
df_sub_p = df2 %>% filter(grepl("^p", label))

# extract info from label
df_label_p = str_split(df_sub_p$label, '_|\\:\\:', simplify=TRUE)
colnames(df_label_p) = c("Category", "rsid", "ref", "alt", "variant", "position", "highAllele", 'genomeCoords')
df_label_p = data.frame(df_label_p)
df_label_p$label = df_sub_p$label

# stack positive and negatives
df_label = smartbind(df_label_p, df_label_n)

# identify unique label matching alt and ref sequences
df_merge = merge(df2, df_label, by="label")
df_merge$label2 = with(df_merge, paste(Category, rsid, ref, alt, position, sep="_"))

# merge alt and ref 
df_match = merge(df_merge %>% 
				filter(variant == 'ref') %>%
				select(Category, rsid, ref, alt, position, label2, DNA1, DNA2, DNA3, RNA1, RNA2, RNA3),
				df_merge %>% 
				filter(variant == 'alt') %>%
				select(label2, DNA1, DNA2, DNA3, RNA1, RNA2, RNA3),
			by = "label2",
			suffixes = c('.ref', '.alt'))

file = paste0(folder_mpra, "mpra_matched_refalt.tsv")
write.table(df_match, file, quote=FALSE, row.names=FALSE)













