# This script reads the yield and persistency gene lists and uses topGO
# to perform an enrichment analysis. Some pointers/useage examples came
# from: https://www.biostars.org/p/9641/. There are also some helpful
# tips in a Slideshare talk by Gomez:
# http://www.slideshare.net/aubombarely/gotermsanalysiswithr.
#
# This is the version of the program for dealing with the intersection
# datasets in which SNP common to all breeds are assembled.

# Load the libraries we need for the analysis.
library(biomaRt)
library(org.Bt.eg.db)
library(topGO)

run_GO_analysis <- function(inputfile, breed1, breed2, breed3, trait, debug=FALSE){

	# Testing...
	if (debug){
		cat("\nParameters used:\n")
		cat(paste("          inputfile:", inputfile, "\n"))
		cat(paste("          breed1   :", breed1, "\n"))
		cat(paste("          breed2   :", breed2, "\n"))
		cat(paste("          breed3   :", breed3, "\n"))
		cat(paste("          trait    :", trait, "\n"))
		cat("\n")
	}

        # Read the data from the input file. The contents should
        # look like this:
        # ENSBTAG00000015100 ATP6V0E1 BovineHD2000001483 0.0049334729
        # ENSBTAG00000004297 ACOXL ARS-BFGL-NGS-38779 0.0046484432
        # ENSBTAG00000009308 TDRD3 BTB-01981334 0.0032824175
        #
	gene_list <- read.table(inputfile, sep=" ", col.names=c("gene_name", "gene_symbol", "snp_name", "breed1_effect", "breed2_effect", "breed3_effect"))
        gene_list$effect <- do.call(pmax, gene_list[4:6])
        cutoff <- min(gene_list$effect)
        print(head(gene_list))

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
		description = paste("topGO analysis of", breed1, breed2, breed3, trait),
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
		fn.prefix = paste(breed1, breed2, breed3, trait, sep="_"),
		useInfo = "all",
		pdfSW = TRUE
	)

	return(go_results)

} # End of function run_GO_analysis().

###
### The run_GO_analysis() function in this file handles ONLY three-breed interactions.
###

# jcole@jcole ~/Documents/AIPL/Genomics/Persistency/networks $ wc -l all_snp_??_??_??_[FMPS]*_shared.txt
#       0 all_snp_ay_bs_ho_Fat_pers_shared.txt
#      39 all_snp_ay_bs_ho_Fat_shared.txt
#      25 all_snp_ay_bs_ho_Milk_pers_shared.txt
#      46 all_snp_ay_bs_ho_Milk_shared.txt
#      11 all_snp_ay_bs_ho_Protein_pers_shared.txt
#      36 all_snp_ay_bs_ho_Protein_shared.txt
#      36 all_snp_ay_bs_ho_SCS_pers_shared.txt
#      40 all_snp_ay_bs_ho_SCS_shared.txt
#      17 all_snp_ay_bs_je_Fat_pers_shared.txt
#      68 all_snp_ay_bs_je_Fat_shared.txt
#      11 all_snp_ay_bs_je_Milk_pers_shared.txt
#      63 all_snp_ay_bs_je_Milk_shared.txt
#      14 all_snp_ay_bs_je_Protein_pers_shared.txt
#      66 all_snp_ay_bs_je_Protein_shared.txt
#       3 all_snp_ay_bs_je_SCS_pers_shared.txt
#      65 all_snp_ay_bs_je_SCS_shared.txt
#      15 all_snp_ay_ho_je_Fat_pers_shared.txt
#      25 all_snp_ay_ho_je_Fat_shared.txt
#      56 all_snp_ay_ho_je_Milk_pers_shared.txt
#      20 all_snp_ay_ho_je_Milk_shared.txt
#       7 all_snp_ay_ho_je_Protein_pers_shared.txt
#      27 all_snp_ay_ho_je_Protein_shared.txt
#       7 all_snp_ay_ho_je_SCS_pers_shared.txt
#      32 all_snp_ay_ho_je_SCS_shared.txt
#      17 all_snp_bs_ho_je_Fat_pers_shared.txt
#      26 all_snp_bs_ho_je_Fat_shared.txt
#      13 all_snp_bs_ho_je_Milk_pers_shared.txt
#      19 all_snp_bs_ho_je_Milk_shared.txt
#      19 all_snp_bs_ho_je_Protein_pers_shared.txt
#      22 all_snp_bs_ho_je_Protein_shared.txt
#       2 all_snp_bs_ho_je_SCS_pers_shared.txt
#      28 all_snp_bs_ho_je_SCS_shared.txt
#     875 total

# Milk yield data
ay_bs_ho_milk_yld <- run_GO_analysis("all_snp_ay_bs_ho_Milk_shared.txt", "ay", "bs", "ho", "Milk")
print(ay_bs_ho_milk_yld)

ay_bs_je_milk_yld <- run_GO_analysis("all_snp_ay_bs_je_Milk_shared.txt", "ay", "bs", "je", "Milk")
print(ay_bs_je_milk_yld)

ay_ho_je_milk_yld <- run_GO_analysis("all_snp_ay_ho_je_Milk_shared.txt", "ay", "ho", "je", "Milk")
print(ay_ho_je_milk_yld)

bs_ho_je_milk_yld <- run_GO_analysis("all_snp_bs_ho_je_Milk_shared.txt", "bs", "ho", "je", "Milk")
print(bs_ho_je_milk_yld)

# Milk pers data
ay_bs_ho_milk_pers <- run_GO_analysis("all_snp_ay_bs_ho_Milk_pers_shared.txt", "ay", "bs", "ho", "Milk_pers")
print(ay_bs_ho_milk_pers)

ay_bs_je_milk_pers <- run_GO_analysis("all_snp_ay_bs_je_Milk_pers_shared.txt", "ay", "bs", "je", "Milk_pers")
print(ay_bs_je_milk_pers)

ay_ho_je_milk_pers <- run_GO_analysis("all_snp_ay_ho_je_Milk_pers_shared.txt", "ay", "ho", "je", "Milk_pers")
print(ay_ho_je_milk_pers)

bs_ho_je_milk_pers <- run_GO_analysis("all_snp_bs_ho_je_Milk_pers_shared.txt", "bs", "ho", "je", "Milk_pers")
print(bs_ho_je_milk_pers)

# Fat yield data
ay_bs_ho_fat_yld <- run_GO_analysis("all_snp_ay_bs_ho_Fat_shared.txt", "ay", "bs", "ho", "Fat")
print(ay_bs_ho_fat_yld)

ay_bs_je_fat_yld <- run_GO_analysis("all_snp_ay_bs_je_Fat_shared.txt", "ay", "bs", "je", "Fat")
print(ay_bs_je_fat_yld)

ay_ho_je_fat_yld <- run_GO_analysis("all_snp_ay_ho_je_Fat_shared.txt", "ay", "ho", "je", "Fat")
print(ay_ho_je_fat_yld)

bs_ho_je_fat_yld <- run_GO_analysis("all_snp_bs_ho_je_Fat_shared.txt", "bs", "ho", "je", "Fat")
print(bs_ho_je_fat_yld)

# Fat pers data
ay_bs_je_fat_pers <- run_GO_analysis("all_snp_ay_bs_je_Fat_pers_shared.txt", "ay", "bs", "je", "Fat_pers")
print(ay_bs_je_fat_pers)
   
ay_ho_je_fat_pers <- run_GO_analysis("all_snp_ay_ho_je_Fat_pers_shared.txt", "ay", "ho", "je", "Fat_pers")
print(ay_ho_je_fat_pers)
   
bs_ho_je_fat_pers <- run_GO_analysis("all_snp_bs_ho_je_Fat_pers_shared.txt", "bs", "ho", "je", "Fat_pers")
print(bs_ho_je_fat_pers)

# Protein yield data
ay_bs_ho_protein_yld <- run_GO_analysis("all_snp_ay_bs_ho_Protein_shared.txt", "ay", "bs", "ho", "Protein")
print(ay_bs_ho_protein_yld)

ay_bs_je_protein_yld <- run_GO_analysis("all_snp_ay_bs_je_Protein_shared.txt", "ay", "bs", "je", "Protein")
print(ay_bs_je_protein_yld)

ay_ho_je_protein_yld <- run_GO_analysis("all_snp_ay_ho_je_Protein_shared.txt", "ay", "ho", "je", "Protein")
print(ay_ho_je_protein_yld)

bs_ho_je_protein_yld <- run_GO_analysis("all_snp_bs_ho_je_Protein_shared.txt", "bs", "ho", "je", "Protein")
print(bs_ho_je_protein_yld)

# Protein pers data
ay_bs_ho_protein_pers <- run_GO_analysis("all_snp_ay_bs_ho_Protein_pers_shared.txt", "ay", "bs", "ho", "Protein_pers")
print(ay_bs_ho_protein_pers)
   
ay_bs_je_protein_pers <- run_GO_analysis("all_snp_ay_bs_je_Protein_pers_shared.txt", "ay", "bs", "je", "Protein_pers")
print(ay_bs_je_protein_pers)
   
ay_ho_je_protein_pers <- run_GO_analysis("all_snp_ay_ho_je_Protein_pers_shared.txt", "ay", "ho", "je", "Protein_pers")
print(ay_ho_je_protein_pers)
   
bs_ho_je_protein_pers <- run_GO_analysis("all_snp_bs_ho_je_Protein_pers_shared.txt", "bs", "ho", "je", "Protein_pers")
print(bs_ho_je_protein_pers)

# SCS yield data
ay_bs_ho_scs_yld <- run_GO_analysis("all_snp_ay_bs_ho_SCS_shared.txt", "ay", "bs", "ho", "SCS")
print(ay_bs_ho_scs_yld)

ay_bs_je_scs_yld <- run_GO_analysis("all_snp_ay_bs_je_SCS_shared.txt", "ay", "bs", "je", "SCS")
print(ay_bs_je_scs_yld)

ay_ho_je_scs_yld <- run_GO_analysis("all_snp_ay_ho_je_SCS_shared.txt", "ay", "ho", "je", "SCS")
print(ay_ho_je_scs_yld)

bs_ho_je_scs_yld <- run_GO_analysis("all_snp_bs_ho_je_SCS_shared.txt", "bs", "ho", "je", "SCS")
print(bs_ho_je_scs_yld)

# SCS pers data
ay_bs_ho_scs_pers <- run_GO_analysis("all_snp_ay_bs_ho_SCS_pers_shared.txt", "ay", "bs", "ho", "SCS_pers")
print(ay_bs_ho_scs_pers)
   
ay_bs_je_scs_pers <- run_GO_analysis("all_snp_ay_bs_je_SCS_pers_shared.txt", "ay", "bs", "je", "SCS_pers")
print(ay_bs_je_scs_pers)
   
ay_ho_je_scs_pers <- run_GO_analysis("all_snp_ay_ho_je_SCS_pers_shared.txt", "ay", "ho", "je", "SCS_pers")
print(ay_ho_je_scs_pers)
   
bs_ho_je_scs_pers <- run_GO_analysis("all_snp_bs_ho_je_SCS_pers_shared.txt", "bs", "ho", "je", "SCS_pers")
print(bs_ho_je_scs_pers)
