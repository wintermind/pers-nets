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

run_GO_analysis <- function(inputfile, breed1, breed2, trait, group, cutoff=0.01, debug=FALSE){

	# Testing...
	if (debug){
		cat("\nParameters used:\n")
		cat(paste("          inputfile:", inputfile, "\n"))
		cat(paste("          breed1   :", breed1, "\n"))
		cat(paste("          breed2   :", breed2, "\n"))
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

	# We have more than one SNP in many genes, so we want to
	# remove duplicates and keep only the entries with the
	# largest effect on each gene. Duplicates cause THE TROUBLE
	# for topGO.
	
	# Sort by gene name and reverse of value
	gene_list_dd <- gene_list[order(gene_list$gene_name, -gene_list$yld_effect),]
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
		description = paste("topGO analysis of", breed, trait),
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
		fn.prefix = paste(breed, trait, sep="_"),
		useInfo = "all",
		pdfSW = TRUE
	)

	return(go_results)

} # End of function run_GO_analysis().

# The cutoffs in the following function calls are based on the allele substitution effect
# for the 200th-largest SNP from each analysis. This is rather arbitrary, so you may change
# it if you like.

# Ayrshire data 
#ay_milk_hp_hy_results <- run_GO_analysis("all_snp_ay_Milk_hi_pers_hi_yld_shared.txt", "ay", "Milk_hp_hy", "yld", 0.7494031725)
#print(ay_milk_hp_hy_results)

#ay_milk_hp_ly_results <- run_GO_analysis("all_snp_ay_Milk_hi_pers_lo_yld_shared.txt", "ay", "    Milk_hp_ly", "pers", 0.0002722999)
#print(ay_milk_hp_ly_results)

ay_milk_lp_hy_results <- run_GO_analysis("all_snp_ay_Milk_lo_pers_hi_yld_shared.txt", "ay", "    Milk_lp_hy", "yld", 0.7494031725)
print(ay_milk_lp_hy_results)
