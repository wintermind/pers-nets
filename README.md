Persistency gene network analysis
=================================
This document describes the computing steps used in a study of genes and gene networks
affecting lactation persistency in four breeds of dairy cattle.

Steps
=====
1. Make a file of the closest gene within 25,000 bp of each SNP

 I modified Derek's original one-liner [sic] to include the gene symbol from his Ensembl 79 file
 and to omit the sex chromosomes.

 `~/bin/closestBed -a chromsome_data_HO.bed -b coordinates/ensembl_gene_annotation.bed -d -t first | perl -lane 'if($F[10] < 25000){print $F[0    ], " ", $F[1], " ", $F[3], " ", $F[7], " ", $F[8];}' > HO_snp_within_25kb_april.txt`

 `jcole@jcole ~/Documents/AIPL/Genomics/Persistency/networks $ wc -l HO_snp_within_25kb_april.txt
   36758 HO_snp_within_25kb_april.txt`

 ```bash
 jcole@jcole ~/Documents/AIPL/Genomics/Persistency/networks $ head HO_snp_within_25kb_april.txt
 chr1 120183 BovineHD0100000035 ENSBTAG00000001753 protein_coding
 chr1 135098 Hapmap43437-BTA-101873 ENSBTAG00000001753 protein_coding
 chr1 158820 BovineHD0100000048 ENSBTAG00000001753 protein_coding
 chr1 183040 BovineHD0100000057 ENSBTAG00000001753 protein_coding
 chr1 208728 BovineHD0100000064 ENSBTAG00000001753 protein_coding
 chr1 267940 ARS-BFGL-NGS-16466 ENSBTAG00000005540 CBX3
 chr1 278952 BovineHD0100000079 ENSBTAG00000005540 CBX3
 chr1 393248 Hapmap34944-BES1_Contig627_1906 ENSBTAG00000001150 KCNE1
 chr1 413044 BovineHD0100000136 ENSBTAG00000001150 KCNE1
 chr1 458718 BovineHD0100000155 ENSBTAG00000013841 C21orf140
 ```

 This file should be useable for all breeds (AY, BS, HO, and JE) because the same SNP
 set is used for all four of them.

2. Now we need to merge the gene names into the breed-trait group datasets. This will
 produce larger files than are strictly necessary, but I'm trading disc space for
 convenience.

 `jcole-mac:networks jcole$ python merge_genes_effects.py`

 The resulting files have four columns, Ensembl gene name, gene symbol, Illumina SNP ID, and
 effect size (in additive genetic SD) for autosomal SNP. I AM NOT INCLUDING X AND PAR
 SNP IN THIS ANALYSIS.

 ```bash
 jcole@jcole ~/Documents/AIPL/Genomics/Persistency/networks $ head all_snp_ho_Milk_april.txt
 ENSBTAG00000026356 DGAT1 ARS-BFGL-NGS-4939 24.026196269
 ENSBTAG00000004761 FOXH1 ARS-BFGL-NGS-57820 21.90689127
 ENSBTAG00000011922 PLEC ARS-BFGL-NGS-107379 15.654802535
 ```

3. Before we can let R do its thing, we need to know what value to use for a cutoff when
deciding which SNP to include in the analysis. I have a bash script that will show the
200th-largest SNP from each dataset.

 ```bash
 jcole@jcole ~/Documents/AIPL/Genomics/Persistency/networks $ ./show_line_200.sh 
 all_snp_ay_Fat_april.txt
 ENSBTAG00000007012 ZNF3 Hapmap39666-BTA-60132 0.0298400073
 all_snp_ay_Fat_pers_april.txt
 ENSBTAG00000000000001623 pseudogene ARS-BFGL-NGS-101744 0.0002583732
 all_snp_ay_Milk_april.txt
 ENSBTAG00000046611 protein_coding Hapmap53602-rs29011081 0.7744312987
 ...
 ```

You can add those values to the runi\_topgo.r script that we're going to use in step 4.

4. The most important thing, probably, is to make sure that you're using the correct
threshold for the cutoff you want to use. Since I'm arbitrarily working with the top
200 effects I need to make sure that run\_topgo.r includes the correct effect size
for each trait-breed combination.

 I like to redirect the output, too:

 `> Rscript run_topgo.r > run_topgo.r.out 2>&1`

 * run\_topgo.r ALSO writes out files to be used with DAVID [http://david.abcc.ncifcrf.gov/]: * 

 ```bash
  >jcole@jcole ~/Documents/AIPL/Genomics/Persistency/networks $ ls *gene*.txt
 ... ay_Fat_pers_genes.txt   ay_Milk_yld_genes.txt      ay_SCS_pers_genes.txt
 ... ay_Fat_yld_genes.txt    ay_Protein_pers_genes.txt  ay_SCS_yld_genes.txt
 ... ay_Milk_pers_genes.txt  ay_Protein_yld_genes.txt   ay_gene_universe.txt
 ```

 * The "\_gene\_universe.txt" file has the complete list of annotated genes within 25 kb *
 * of the SNP used in the analysis, and the other files have the SNP with the 200 *
 * largest effects for each trait (the gene set). *

5. Identify the SNP that have different effects on yield and persistency:

   1. Large effects on both persistency and yield
   2. Small effects in yield and large effects on persistency
   3. Large effects on yield and small effects on persistency

 This is done with the find\_intersection.sh script:

 ```bash
 > bash sh find_intersection.sh > find_intersection.out & tail -f find_intersection.out
 ... Breed: ay
 ... 	Trait: Milk
 ... 		Preparing filenames.
 ... 		Large pers, large yld
 ... 			Sorting input files.
 ... 			Computing intersections.
 ... <etc.>
 ```

 Now we're going to use a version of the run\_topgo.r script to work on the intersections. We're going to
 go ahead and use the cutoff value for (1) the yld trait in the high-high analysis, (2) the yld trait in
 the low pers-high yld analysis, and (3) the pers trait in the high pers-low yld analysis.

 The run\_topgo\_intersections.r script will automatically assign the largest effect size (yield or persistency)
 to a gene-SNP pair and use that in the analysis. It also retains only ONE record for each gene, the one
 associated with the largest effect size, because topGO cannot handle a case where there is more than one
 SNP per gene.

* You have to look in find\_intersection.out to see what analyses produced empty sets (no genes in the *
* intersection) so that you can update run\_topgo\_intersections.r accordingly. If you don't, R will *
* throw errors on the empty files. *

6. Now we're going to look at SNP from among the 500 largest effects in each breed that are in common
to two or more breeds (the intersections). We're using the top 500 instead of the top 200 based on a
purely empirical estimate -- below this and you get few/no intersections, and above this and you get
a lot.

 ```bash 
 > bash find_intersection_breeds.sh
 ... Computing intersections of the top 500 SNP for all breed combinations.
 ... 	Trait: Milk
 ... 		Sorting and subsetting top 500 SNP 
 ... 		(Note: There may be more than 1	SNP per gene)
 ... 		AY-BS
 ... 			56 all_snp_ay_bs_Milk_shared.txt genes
 ... 			48 all_snp_ay_bs_Milk_pers_shared.txt genes
 ... <etc.>
 ```

* You __MUST__ use bash to execute this script, not sh!!! *

7. Now I have to figure out what all of this means...
