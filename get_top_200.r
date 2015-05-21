# Read the data from the input file. The contents should
# look like this:
# ENSBTAG00000015100 ATP6V0E1 BovineHD2000001483 0.0049334729
# ENSBTAG00000004297 ACOXL ARS-BFGL-NGS-38779 0.0046484432
# ENSBTAG00000009308 TDRD3 BTB-01981334 0.0032824175
#
gene_list <- read.table("all_snp_ay_Milk_april.txt", sep=" ", col.names=c("gene_name", "gene_symbol", "snp_name", "effect"))
print(head(gene_list))

# We have more than one SNP in many genes, so we want to
# remove duplicates and keep only the entries with the
# largest effect on each gene. Duplicates cause THE TROUBLE
# for topGO.
 
# Sort by gene name and reverse of value
gene_list_dd <- gene_list[order(gene_list$gene_name, -gene_list$effect),]

# Write JUST the Ensembl IDs to a file to create the universe for DAVID.
write.table(gene_list_dd$gene_name, 'ay_gene_universe.txt', row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")

# We want to pull out the genes associated with the 200 largest SNP effects
# to put them into DAVID.
write.table(gene_list_dd[1:200,]$gene_name, 'ay_Milk_genes.txt', row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")

# Take the first row for each gene
#gene_list_u <- gene_list_dd[ !duplicated(gene_list_dd$gene_name), ]
