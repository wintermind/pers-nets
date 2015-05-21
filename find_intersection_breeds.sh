#!/bin/bash
# Create a new text file that contains the Ensembl gene names for all
# genes that are associated with the ${topn} largest SNP effects for each
# pair of traits (yield and persistency).
#
# How do the mysterious join statements work? See here:
# http://stackoverflow.com/questions/2619562/joining-multiple-fields-using-unix-join.

# aybs, ayho, ayje, bsho, bsje, hoje
# aybsho, aybsje, ayhoje, bshoje
# aybshoje

traits="Milk Fat Protein SCS"
breeds="ay bs ho je"
topn=500

printf "\nComputing intersections of the top ${topn} SNP for all breed combinations."
for trait in $traits; do
	printf "\n\tTrait: ${trait}"
	## First, we're going to do the sorting and subsetting common to all
	## sets so that we don't keep repeating operations that only need to
	## be done once.
	printf "\n\t\tSorting and subsetting top ${topn} SNP \n\t\t(Note: There may be more than 1 SNP per gene)"
	for breed in $breeds; do
		# Do the sorting
		sort -g -r -k 3 "all_snp_${breed}_${trait}_april.txt" -o "all_snp_${breed}_${trait}.txt"
		sort -g -r -k 3 "all_snp_${breed}_${trait}_pers_april.txt" -o "all_snp_${breed}_${trait}_pers.txt"
		# Keep the top ${topn} from each set
		head -n ${topn} "all_snp_${breed}_${trait}_april.txt" > "all_snp_${breed}_${trait}.txt.${topn}"
		head -n ${topn} "all_snp_${breed}_${trait}_pers_april.txt" > "all_snp_${breed}_${trait}_pers.txt.${topn}"
		# Sort the subsets by gene ID
		sort -k 1 "all_snp_${breed}_${trait}.txt.${topn}" -o "all_snp_${breed}_${trait}.txt.${topn}"
		sort -k 1 "all_snp_${breed}_${trait}_pers.txt.${topn}" -o "all_snp_${breed}_${trait}_pers.txt.${topn}"
	done

	##
	## Pairwise intersections
	##

	## AY-BS
	printf "\n\t\tAY-BS"
	join -j1 -o 1.2,1.3,1.4,2.4 <(<"all_snp_ay_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_bs_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_bs_${trait}_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_bs_${trait}_shared.txt") genes"
	join -j1 -o 1.2,1.3,1.4,2.4 <(<"all_snp_ay_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_bs_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_bs_${trait}_pers_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_bs_${trait}_pers_shared.txt") genes"

	## AY-HO
	printf "\n\t\tAY-HO"
	join -j1 -o 1.2,1.3,1.4,2.4 <(<"all_snp_ay_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_ho_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_ho_${trait}_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_ho_${trait}_shared.txt") genes"
	join -j1 -o 1.2,1.3,1.4,2.4 <(<"all_snp_ay_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_ho_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_ho_${trait}_pers_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_ho_${trait}_pers_shared.txt") genes"

	## AY-JE
	printf "\n\t\tAY-JE"
	join -j1 -o 1.2,1.3,1.4,2.4 <(<"all_snp_ay_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_je_${trait}_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_je_${trait}_shared.txt") genes"
	join -j1 -o 1.2,1.3,1.4,2.4 <(<"all_snp_ay_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_je_${trait}_pers_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_je_${trait}_pers_shared.txt") genes"

	## BS-HO
	printf "\n\t\tBS-HO"
	join -j1 -o 1.2,1.3,1.4,2.4 <(<"all_snp_bs_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_ho_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_bs_ho_${trait}_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_bs_ho_${trait}_shared.txt") genes"
	join -j1 -o 1.2,1.3,1.4,2.4 <(<"all_snp_bs_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_ho_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_bs_ho_${trait}_pers_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_bs_ho_${trait}_pers_shared.txt") genes"

	## BS-JE
	printf "\n\t\tBS-JE"
	join -j1 -o 1.2,1.3,1.4,2.4 <(<"all_snp_bs_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_bs_je_${trait}_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_bs_je_${trait}_shared.txt") genes"
	join -j1 -o 1.2,1.3,1.4,2.4 <(<"all_snp_bs_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_bs_je_${trait}_pers_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_bs_je_${trait}_pers_shared.txt") genes"

	## HO-JE
    printf "\n\t\tHO-JE"
	join -j1 -o 1.2,1.3,1.4,2.4 <(<"all_snp_ho_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ho_je_${trait}_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ho_je_${trait}_shared.txt") genes"
	join -j1 -o 1.2,1.3,1.4,2.4 <(<"all_snp_ho_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ho_je_${trait}_pers_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ho_je_${trait}_pers_shared.txt") genes"

	##
	## Three-way interactions
	##

	## AY-BS-HO
	printf "\n\t\tAY-BS-HO"
	join -j1 -o 1.2,1.3,1.4,1.5,2.4 <(<"all_snp_ay_bs_${trait}_shared.txt" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_ho_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_bs_ho_${trait}_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_bs_ho_${trait}_shared.txt") genes"
	join -j1 -o 1.2,1.3,1.4,1.5,2.4 <(<"all_snp_ay_bs_${trait}_pers_shared.txt" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_ho_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_bs_ho_${trait}_pers_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_bs_ho_${trait}_pers_shared.txt") genes"

	## AY-BS-JE
	printf "\n\t\tAY-BS-JE"
	join -j1 -o 1.2,1.3,1.4,1.5,2.4 <(<"all_snp_ay_bs_${trait}_shared.txt" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_bs_je_${trait}_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_bs_je_${trait}_shared.txt") genes"
	join -j1 -o 1.2,1.3,1.4,1.5,2.4 <(<"all_snp_ay_bs_${trait}_pers_shared.txt" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_bs_je_${trait}_pers_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_bs_je_${trait}_pers_shared.txt") genes"

	## AY-HO-JE
	printf "\n\t\tAY-HO-JE"
	join -j1 -o 1.2,1.3,1.4,1.5,2.4 <(<"all_snp_ay_ho_${trait}_shared.txt" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_ho_je_${trait}_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_ho_je_${trait}_shared.txt") genes"
	join -j1 -o 1.2,1.3,1.4,1.5,2.4 <(<"all_snp_ay_ho_${trait}_pers_shared.txt" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_ho_je_${trait}_pers_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_ho_je_${trait}_pers_shared.txt") genes"

	## BS-HO-JE
	printf "\n\t\tBS-HO-JE"
	join -j1 -o 1.2,1.3,1.4,1.5,2.4 <(<"all_snp_bs_ho_${trait}_shared.txt" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_bs_ho_je_${trait}_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_bs_ho_je_${trait}_shared.txt") genes"
	join -j1 -o 1.2,1.3,1.4,1.5,2.4 <(<"all_snp_bs_ho_${trait}_pers_shared.txt" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_bs_ho_je_${trait}_pers_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_bs_ho_je_${trait}_pers_shared.txt") genes"

	##
	## Four-way interactions
    ##

	## AY-BS-HO-JE
	printf "\n\t\tAY-BS-HO-JE"
    join -j1 -o 1.2,1.3,1.4,1.5,1.6,2.4 <(<"all_snp_ay_bs_ho_${trait}_shared.txt" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_bs_ho_je_${trait}_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_bs_ho_je_${trait}_shared.txt") genes"
    join -j1 -o 1.2,1.3,1.4,1.5,1.6,2.4 <(<"all_snp_ay_bs_ho_${trait}_pers_shared.txt" awk '{print $1"-"$2" "$0}' | sort -k1,1) <(<"all_snp_je_${trait}_pers.txt.${topn}" awk '{print $1"-"$2" "$0}' | sort -k1,1) > "all_snp_ay_bs_ho_je_${trait}_pers_shared.txt"
	printf "\n\t\t\t$(wc -l "all_snp_ay_bs_ho_je_${trait}_pers_shared.txt") genes"

done
printf "\n"
