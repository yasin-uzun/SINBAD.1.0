

#Read regions
print('Reading gene coordinates')
df_gene_annot = read_region_annot(gene_annot_file, format_file)
head(df_gene_annot)

df_gene_annot_2kb = df_gene_annot
df_gene_annot_2kb$start = df_gene_annot$start - 2000
df_gene_annot_2kb$start[df_gene_annot_2kb$start < 0 ] = 0
df_gene_annot_2kb$end = df_gene_annot$end + 2000

head(df_gene_annot_2kb)
genes_2kb_annot_file = "/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/dependencies//annot//hg38/regions.all_genes_plus_2Kb.bed"
write.table(df_gene_annot_2kb, genes_2kb_annot_file, 
            row.names = F, col.names = F, quote = F, sep = '\t')




print('Computing promoters')

head(df_gene_annot)
annot_promoters = df_gene_annot
head(annot_promoters)
attach(df_gene_annot)
annot_promoters$start[strand == '+'] = df_gene_annot$start[strand == '+'] - 2000
annot_promoters$end[strand == '+'] = df_gene_annot$start[strand == '+'] + 500

annot_promoters$end[strand == '-'] = df_gene_annot$end[strand == '-'] + 2000
annot_promoters$start[strand == '-'] = df_gene_annot$end[strand == '-'] - 500
detach(df_gene_annot)
head(annot_promoters)

annot_promoters$start[annot_promoters$start < 0] = 0

head(annot_promoters)

promoters_file = "/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/dependencies//annot//hg38/regions.TSS_up2Kb_down500bp.bed"
write.table(annot_promoters, promoters_file, 
            row.names = F, col.names = F, quote = F, sep = '\t')





head(df_gene_annot)
annot_promoters = df_gene_annot
head(annot_promoters)
attach(df_gene_annot)
annot_promoters$start[strand == '+'] = df_gene_annot$start[strand == '+'] - 2000
annot_promoters$end[strand == '+'] = df_gene_annot$start[strand == '+'] + 2000

annot_promoters$end[strand == '-'] = df_gene_annot$end[strand == '-'] + 2000
annot_promoters$start[strand == '-'] = df_gene_annot$end[strand == '-'] - 2000
detach(df_gene_annot)
head(annot_promoters)

annot_promoters$start[annot_promoters$start < 0] = 0

head(annot_promoters)

promoters_file = "/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/dependencies//annot//hg38/regions.TSS_up2Kb_down2Kb.bed"
write.table(annot_promoters, promoters_file, 
            row.names = F, col.names = F, quote = F, sep = '\t')








