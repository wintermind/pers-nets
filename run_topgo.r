# This script reads the yield and persistency gene lists and uses topGO
# to perform an enrichment analysis. Some pointers/useage examples came
# from: https://www.biostars.org/p/9641/. There are also some helpful
# tips in a Slideshare talk by Gomez:
# http://www.slideshare.net/aubombarely/gotermsanalysiswithr.

# Load the libraries we need for the analysis.
library(biomaRt)
library(org.Bt.eg.db)
library(topGO)

run_GO_analysis <- function(inputfile, breed, trait, group, cutoff=0.01, debug=FALSE){

	# Testing...
	if (debug){
		cat("\nParameters used:\n")
		cat(paste("          inputfile:", inputfile, "\n"))
		cat(paste("          breed    :", breed, "\n"))
		cat(paste("          trait    :", trait, "\n"))
        	cat(paste("          group    :", group, "\n"))
		cat(paste("          cutoff   :", cutoff, "\n"))
		cat("\n")
	}

	# Read the data from the input file. The contents should
	# look like this:
	# ENSBTAG00000015100 ATP6V0E1 BovineHD2000001483 0.0049334729
	# ENSBTAG00000004297 ACOXL ARS-BFGL-NGS-38779 0.0046484432
	# ENSBTAG00000009308 TDRD3 BTB-01981334 0.0032824175
	#
	gene_list <- read.table(inputfile, sep=" ", col.names=c("gene_name", "gene_symbol", "snp_name", "effect"))
	print(head(gene_list))

	# Create the gene list files we need for DAVID.
	write_david_files(inputfile, breed, trait, group, debug)

	# We have more than one SNP in many genes, so we want to
	# remove duplicates and keep only the entries with the
	# largest effect on each gene. Duplicates cause THE TROUBLE
	# for topGO.
	
	# Sort by gene name and reverse of value
	gene_list_dd <- gene_list[order(gene_list$gene_name, -gene_list$effect),]
	# Take the first row for each gene
	gene_list_u <- gene_list_dd[ !duplicated(gene_list_dd$gene_name), ]

	# Now we need to pull out the vector of effects and name it so that
	# topGO won't complain.
	all_genes <- t(gene_list_u$effect)
	names(all_genes) <- t(gene_list_u$gene_name)
	cat("\n")
	cat(paste("        rows in all_genes :", nrow(all_genes), "\n"))
	cat(paste("        cols in all_genes :", ncol(all_genes), "\n"))

	# Define the function for selecting interesting genes from the
	# gene universe.
	geneSelFunc <- function(score) {
		# We're going to read (not change!) the variable "cutoff" from
		# the global namespace. Don't tell anyone on Stack Overflow...
		# The reason I did this is that the cutoffs may vary by breed
		# and trait.
		#cat(paste("\n[geneSelFun]: cutoff = ", cutoff,"\n"))
		return(score >= cutoff)
	}

	sel_genes <- geneSelFunc(all_genes)
	names(sel_genes) <- t(gene_list_u$gene_name)
	print(sel_genes[1:10])

	# Now to run topGO. Note that the geneSel argument MUST be changed if you
	# provide a difference value for "effect" in the input file. I'm using
	# the magnitude of the allele substitution effect expressed in additive
	# genetic standard deviations, but you could use p-values or -log10(p-val)
	# instead.
	GOdata <- new("topGOdata",
		ontology = "BP",
		allGenes = all_genes,
		geneSel = geneSelFunc,
		description = paste("topGO analysis of", breed, trait, group),
		annot = annFUN.org,
		mapping="org.Bt.eg.db",
		ID="Ensembl"
	)

	# Perform Fisher's exact test 
	rFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

	# Compute a test statistic for the Fisher's exact test results
	test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
	resultWeight <- getSigGroups(GOdata, test.stat)

	# Perform a Kolmogorov-Smirnov test
	rKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")

	# Construct the table of results.
	go_results <- GenTable(GOdata,
		classic = rFisher,
		KS = rKS,
		weight = resultWeight,
		orderBy = "weight",
		ranksOf = "classic",
		topNodes = 20
	)

	# The subgraph induced by the 5 most significant GO terms as identified by the elim
	# algorithm. Significant nodes are represented as rectangles.
	showSigOfNodes(GOdata, score(rKS), firstSigNodes = 5, useInfo = 'all')
	printGraph(GOdata,
		rKS,
		firstSigNodes = 5,
		fn.prefix = paste(breed, trait, group, sep="_"),
		useInfo = "all",
		pdfSW = TRUE
	)

	return(go_results)

} # End of function run_GO_analysis().

write_david_files <- function(inputfile, breed, trait, group, debug=FALSE){

	# Testing...
	if (debug){
		cat("\nParameters used:\n")
		cat(paste("          inputfile:", inputfile, "\n"))
		cat(paste("          breed    :", breed, "\n"))
		cat(paste("          trait    :", trait, "\n"))
		cat(paste("          group    :", group, "\n"))
		cat("\n")
	}

	# Read the data from the input file. The contents should
	# look like this:
	# ENSBTAG00000015100 ATP6V0E1 BovineHD2000001483 0.0049334729
	# ENSBTAG00000004297 ACOXL ARS-BFGL-NGS-38779 0.0046484432
	# ENSBTAG00000009308 TDRD3 BTB-01981334 0.0032824175
	#
	gene_list <- read.table("all_snp_ay_Milk_april.txt", sep=" ", col.names=c("gene_name", "gene_symbol", "snp_name", "effect"))

	universe_file <- paste(breed,"gene_universe.txt",sep="_")
	david_file <- paste(breed,trait,group,"genes.txt",sep="_")

	# We have more than one SNP in many genes, so we want to
	# remove duplicates and keep only the entries with the
	# largest effect on each gene. Duplicates cause THE TROUBLE
	# for topGO.
 
	# Sort by gene name and reverse of value
	gene_list_dd <- gene_list[order(gene_list$gene_name, -gene_list$effect),]

	# Write JUST the Ensembl IDs to a file to create the universe for DAVID.
	write.table(gene_list_dd$gene_name, universe_file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
 
	# We want to pull out the genes associated with the 200 largest SNP effects
	# to put them into DAVID.
	write.table(gene_list_dd[1:200,]$gene_name, david_file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")

} # End of function write_david_files().

# The cutoffs in the following function calls are based on the allele substitution effect
# for the 200th-largest SNP from each analysis. This is rather arbitrary, so you may change
# it if you like.

# Ayrshire data 
ay_milk_results <- run_GO_analysis("all_snp_ay_Milk_april.txt", "ay", "Milk", "yld", 0.7744312987)
print(ay_milk_results)
ay_milk_pers_results <- run_GO_analysis("all_snp_ay_Milk_pers_april.txt", "ay", "Milk", "pers", 0.0002815349)
print(ay_milk_pers_results)
ay_fat_results <- run_GO_analysis("all_snp_ay_Fat_april.txt", "ay", "Fat", "yld", 0.0298400073)
print(ay_fat_results)
ay_fat_pers_results <- run_GO_analysis("all_snp_ay_Fat_pers_april.txt", "ay", "Fat", "pers", 0.0002583732)
print(ay_fat_pers_results)
ay_protein_results <- run_GO_analysis("all_snp_ay_Protein_april.txt", "ay", "Protein", "yld", 0.022452541)
print(ay_protein_results)
ay_protein_pers_results <- run_GO_analysis("all_snp_ay_Protein_pers_april.txt", "ay", "Protein", "pers", 0.0002931191)
print(ay_protein_pers_results)
ay_scs_results <- run_GO_analysis("all_snp_ay_SCS_april.txt", "ay", "SCS", "yld", 0.00028375)
print(ay_scs_results)
ay_scs_pers_results <- run_GO_analysis("all_snp_ay_SCS_pers_april.txt", "ay", "SCS", "pers", 0.0001466711)
print(ay_scs_pers_results)

# Brown Swiss data 
bs_milk_results <- run_GO_analysis("all_snp_bs_Milk_april.txt", "bs", "Milk", "yld", 1.6681821695)
print(bs_milk_results)
bs_milk_pers_results <- run_GO_analysis("all_snp_bs_Milk_pers_april.txt", "bs", "Milk", "pers", 0.0003839432)
print(bs_milk_pers_results)
bs_fat_results <- run_GO_analysis("all_snp_bs_Fat_april.txt", "bs", "Fat", "yld", 0.0606780749)
print(bs_fat_results)
bs_fat_pers_results <- run_GO_analysis("all_snp_bs_Fat_pers_april.txt", "bs", "Fat", "pers", 0.0003495239)
print(bs_fat_pers_results)
bs_protein_results <- run_GO_analysis("all_snp_bs_Protein_april.txt", "bs", "Protein", "yld", 0.0471150877)
print(bs_protein_results)
bs_protein_pers_results <- run_GO_analysis("all_snp_bs_Protein_pers_april.txt", "bs", "Protein", "pers", 0.0003286103)
print(bs_protein_pers_results)
bs_scs_results <- run_GO_analysis("all_snp_bs_SCS_april.txt", "bs", "SCS", "yld", 0.0005218492)
print(bs_scs_results)
bs_scs_pers_results <- run_GO_analysis("all_snp_bs_SCS_pers_april.txt", "bs", "SCS", "pers", 0.0002322627)
print(bs_scs_pers_results)

# Holstein data
ho_milk_results <- run_GO_analysis("all_snp_ho_Milk_april.txt", "ho", "Milk", "yld", 2.9871631466)
print(ho_milk_results)
ho_milk_pers_results <- run_GO_analysis("all_snp_ho_Milk_pers_april.txt", "ho", "Milk", "pers", 0.001868695)
print(ho_milk_pers_results)
ho_fat_results <- run_GO_analysis("all_snp_ho_Fat_april.txt", "ho", "Fat", "yld", 0.1026451294)
print(ho_fat_results)
ho_fat_pers_results <- run_GO_analysis("all_snp_ho_Fat_pers_april.txt", "ho", "Fat", "pers", 0.0015198632)
print(ho_fat_pers_results)
ho_protein_results <- run_GO_analysis("all_snp_ho_Protein_april.txt", "ho", "Protein", "yld", 0.075492098)
print(ho_protein_results)
ho_protein_pers_results <- run_GO_analysis("all_snp_ho_Protein_pers_april.txt", "ho", "Protein", "pers", 0.0014462046)
print(ho_protein_pers_results)
ho_scs_results <- run_GO_analysis("all_snp_HO_SCS_april.txt", "ho", "SCS", "yld", 0.0009768092)
print(ho_scs_results)
ho_scs_pers_results <- run_GO_analysis("all_snp_ho_SCS_pers_april.txt", "ho", "SCS", "pers", 0.0007462432)
print(ho_scs_pers_results)

# Jersey data 
je_milk_results <- run_GO_analysis("all_snp_je_Milk_april.txt", "je", "Milk", "yld", 1.995651183)
print(je_milk_results)
je_milk_pers_results <- run_GO_analysis("all_snp_je_Milk_pers_april.txt", "je", "Milk", "pers", 0.0012678466)
print(je_milk_pers_results)
je_fat_results <- run_GO_analysis("all_snp_je_Fat_april.txt", "je", "Fat", "yld", 0.0765808207)
print(je_fat_results)
je_fat_pers_results <- run_GO_analysis("all_snp_je_Fat_pers_april.txt", "je", "Fat", "pers", 0.0010669253)
print(je_fat_pers_results)
je_protein_results <- run_GO_analysis("all_snp_je_Protein_april.txt", "je", "Protein", "yld", 0.0617138885)
print(je_protein_results)
je_protein_pers_results <- run_GO_analysis("all_snp_je_Protein_pers_april.txt", "je", "Protein", "pers", 0.0011200394)
print(je_protein_pers_results)
je_scs_results <- run_GO_analysis("all_snp_je_SCS_april.txt", "je", "SCS", "yld", 0.0005449717)
print(je_scs_results)
je_scs_pers_results <- run_GO_analysis("all_snp_je_SCS_pers_april.txt", "je", "SCS", "pers", 0.0004937861)
print(je_scs_pers_results)
