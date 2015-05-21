# This script reads the yield and persistency gene lists and uses topGO
# to perform an enrichment analysis. Some pointers/useage examples came
# from: https://www.biostars.org/p/9641/. There are also some helpful
# tips in a Slideshare talk by Gomez:
# http://www.slideshare.net/aubombarely/gotermsanalysiswithr.
#
# This is the version of the program for dealing with the intersection
# datasets.

# Load the libraries we need for the analysis.
library(biomaRt)
library(org.Bt.eg.db)
library(topGO)

run_GO_analysis <- function(inputfile, breed, trait, debug=FALSE){

	# Testing...
	if (debug){
		cat("\nParameters used:\n")
		cat(paste("          inputfile:", inputfile, "\n"))
		cat(paste("          breed    :", breed, "\n"))
		cat(paste("          trait    :", trait, "\n"))
		cat("\n")
	}

        # Read the data from the input file. The contents should
        # look like this:
	# ==> all_snp_ho_hh_Milk_shared.txt <==
	# ENSBTAG00000000084 MRPS30 BovineHD2000008833 2.4296776706 0.0018258655
	# ENSBTAG00000000219 GRIN2B Hapmap51915-BTA-74618 2.6781408512 0.0017316314
	# ENSBTAG00000000395 ADAMTS20 BovineHD0500010672 2.8259553695 0.0019333065
        gene_list <- read.table(inputfile, sep=" ", col.names=c("gene_name", "gene_symbol", "snp_name", "yld_effect", "pers_effect"))
	gene_list$effect <- do.call(pmax, gene_list[4:5])
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

### You need to look in find_intersections.out to see what intersections had no members.
### If you try and run this on empty files, well, R gets cross.

## SNP with large effects on yield and large effects on persistency

# Milk data
ay_hh_milk_results <- run_GO_analysis("all_snp_ay_hh_Milk_shared.txt", "ay", "Milk")
print(ay_hh_milk_results)

bs_hh_milk_results <- run_GO_analysis("all_snp_bs_hh_Milk_shared.txt", "bs", "Milk")
print(bs_hh_milk_results)

ho_hh_milk_results <- run_GO_analysis("all_snp_ho_hh_Milk_shared.txt", "ho", "Milk")
print(ho_hh_milk_results)

je_hh_milk_results <- run_GO_analysis("all_snp_je_hh_Milk_shared.txt", "je", "Milk")
print(je_hh_milk_results)

# Fat data
ay_hh_fat_results <- run_GO_analysis("all_snp_ay_hh_Fat_shared.txt", "ay", "Fat")
print(ay_hh_fat_results)

bs_hh_fat_results <- run_GO_analysis("all_snp_bs_hh_Fat_shared.txt", "bs", "Fat")
print(bs_hh_fat_results)

ho_hh_fat_results <- run_GO_analysis("all_snp_ho_hh_Fat_shared.txt", "ho", "Fat")
print(ho_hh_fat_results)

je_hh_fat_results <- run_GO_analysis("all_snp_je_hh_Fat_shared.txt", "je", "Fat")
print(je_hh_fat_results)

# Protein data
ay_hh_protein_results <- run_GO_analysis("all_snp_ay_hh_Protein_shared.txt", "ay", "Protein")
print(ay_hh_milk_results)

bs_hh_protein_results <- run_GO_analysis("all_snp_bs_hh_Protein_shared.txt", "bs", "Protein")
print(bs_hh_protein_results)

ho_hh_protein_results <- run_GO_analysis("all_snp_ho_hh_Protein_shared.txt", "ho", "Protein")
print(ho_hh_protein_results)

je_hh_protein_results <- run_GO_analysis("all_snp_je_hh_Protein_shared.txt", "je", "Protein")
print(je_hh_protein_results)

# SCS data
ay_hh_scs_results <- run_GO_analysis("all_snp_ay_hh_SCS_shared.txt", "ay", "SCS")
print(ay_hh_scs_results)

bs_hh_scs_results <- run_GO_analysis("all_snp_bs_hh_SCS_shared.txt", "bs", "SCS")
print(bs_hh_scs_results)

ho_hh_scs_results <- run_GO_analysis("all_snp_ho_hh_SCS_shared.txt", "ho", "SCS")
print(ho_hh_scs_results)

je_hh_scs_results <- run_GO_analysis("all_snp_je_hh_SCS_shared.txt", "je", "SCS")
print(je_hh_scs_results)

## SNP with large effects on yield and small effects on persistency

# Milk data
ho_hl_milk_results <- run_GO_analysis("all_snp_ho_hl_Milk_shared.txt", "ho", "Milk")
print(ho_hl_milk_results)

# Fat data
ho_hl_fat_results <- run_GO_analysis("all_snp_ho_hl_Fat_shared.txt", "ho", "Fat")
print(ho_hl_fat_results)

# Protein data
ho_hl_protein_results <- run_GO_analysis("all_snp_ho_hl_Protein_shared.txt", "ho", "Protein")
print(ho_hl_protein_results)

# SCS data
ho_hl_scs_results <- run_GO_analysis("all_snp_ho_hl_SCS_shared.txt", "ho", "SCS")
print(ho_hl_scs_results)

## SNP with small effects on yield and large effects on persistency

# Milk data
ho_lh_milk_results <- run_GO_analysis("all_snp_ho_lh_Milk_shared.txt", "ho", "Milk")
print(ho_lh_milk_results)

# Fat data

# Protein data
ho_lh_protein_results <- run_GO_analysis("all_snp_ho_lh_Protein_shared.txt", "ho", "Protein")
print(ho_lh_protein_results)

je_lh_protein_results <- run_GO_analysis("all_snp_je_lh_Protein_shared.txt", "je", "Protein")
print(je_lh_protein_results)

# SCS data
ho_lh_scs_results <- run_GO_analysis("all_snp_ho_lh_SCS_shared.txt", "ho", "SCS")
print(ho_lh_scs_results)
