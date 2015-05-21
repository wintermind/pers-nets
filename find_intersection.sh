#!/bin/bash
# Create a new text file that contains the Ensembl gene names for all
# genes that are associated with the ${topn} largest SNP effects for each
# pair of traits (yield and persistency).
#
# How do the mysterious join statements work? See here:
# http://stackoverflow.com/questions/2619562/joining-multiple-fields-using-unix-join.

traits="Milk Fat Protein SCS"
breeds="ay bs ho je"
topn=500

printf "        Finding SNP among the top ${topn} with large effects on yield and persistency."
for trait in $traits; do
	printf "\n        Trait: ${trait}"
	## First, we're going to do the sorting and subsetting common to all
	## sets so that we don't keep repeating operations that only need to
	## be done once.
	printf "\n                Sorting and subsetting top ${topn} SNP"
        printf "\n                (Note: There may be more than 1 SNP per gene)"
	for breed in $breeds; do
		# Do the sorting
		sort -g -r -k 4 "all_snp_${breed}_${trait}_april.txt" -o "all_snp_${breed}_${trait}_april.txt"
		sort -g -r -k 4 "all_snp_${breed}_${trait}_pers_april.txt" -o "all_snp_${breed}_${trait}_pers_april.txt"
		# Keep the top ${topn} from each set
		head -n ${topn} "all_snp_${breed}_${trait}_april.txt" > "all_snp_${breed}_${trait}_april.txt.${topn}"
		head -n ${topn} "all_snp_${breed}_${trait}_pers_april.txt" > "all_snp_${breed}_${trait}_pers_april.txt.${topn}"
		# Sort the subsets by gene ID
		sort -k 1 "all_snp_${breed}_${trait}_april.txt.${topn}" -o "all_snp_${breed}_${trait}_april.txt.${topn}"
		sort -k 1 "all_snp_${breed}_${trait}_pers_april.txt.${topn}" -o "all_snp_${breed}_${trait}_pers_april.txt.${topn}"
		# Compute the intersections
		printf "\n                ${breed} -- ${trait}"
		join -j1 -o 1.2,1.3,1.4,1.5,2.5 <(<"all_snp_${breed}_${trait}_april.txt.${topn}" awk '{print $1"-"$2"-"$3" "$0}' | sort -k1,1) <(<"all_snp_${breed}_${trait}_pers_april.txt.${topn}" awk '{print $1"-"$2"-"$3" "$0}' | sort -k1,1) > "all_snp_${breed}_hh_${trait}_shared.txt"
		printf "\n                        $(wc -l "all_snp_${breed}_hh_${trait}_shared.txt") genes"
		# Clean-up
		rm "all_snp_${breed}_${trait}_april.txt.${topn}"
		rm "all_snp_${breed}_${trait}_pers_april.txt.${topn}"
	done
done
printf "\n"

printf "        Finding SNP among the top ${topn} with large effects on yield and small effects on persistency."
for trait in $traits; do
	printf "\n        Trait: ${trait}"
    	## First, we're going to do the sorting and subsetting common to all
     	## sets so that we don't keep repeating operations that only need to
      	## be done once.
       	printf "\n                Sorting and subsetting top ${topn} SNP"
        printf "\n                (Note: There may be more than 1 SNP per gene)"
	for breed in $breeds; do
		# Do the sorting
	        sort -g -r -k 4 "all_snp_${breed}_${trait}_april.txt" -o "all_snp_${breed}_${trait}_april.txt"
	        sort -g -k 4 "all_snp_${breed}_${trait}_pers_april.txt" -o "all_snp_${breed}_${trait}_pers_april.txt"
	        # Keep the top ${topn} from each set
	        head -n ${topn} "all_snp_${breed}_${trait}_april.txt" > "all_snp_${breed}_${trait}_april.txt.${topn}"
	        head -n ${topn} "all_snp_${breed}_${trait}_pers_april.txt" > "all_snp_${breed}_${trait}_pers_april.txt.${topn}"
	        # Sort the subsets by gene ID
	        sort -k 1 "all_snp_${breed}_${trait}_april.txt.${topn}" -o "all_snp_${breed}_${trait}_april.txt.${topn}"
	        sort -k 1 "all_snp_${breed}_${trait}_pers_april.txt.${topn}" -o "all_snp_${breed}_${trait}_pers_april.txt.${topn}"
	        # Compute the intersections
	        printf "\n                ${breed} -- ${trait}"
		join -j1 -o 1.2,1.3,1.4,1.5,2.5 <(<"all_snp_${breed}_${trait}_april.txt.${topn}" awk '{print $1"-"$2"-"$3" "$0}' | sort -k1,1) <(<"all_snp_${breed}_${trait}_pers_april.txt.${topn}" awk '{print $1"-"$2"-"$3" "$0}' | sort -k1,1) > "all_snp_${breed}_hl_${trait}_shared.txt"
		printf "\n                        $(wc -l "all_snp_${breed}_hl_${trait}_shared.txt") genes"
		# Clean-up
		rm "all_snp_${breed}_${trait}_april.txt.${topn}"
		rm "all_snp_${breed}_${trait}_pers_april.txt.${topn}"
	done
done
printf "\n"

 printf "        Finding SNP among the top ${topn} with small effects on yield and large effects on persistency."
 for trait in $traits; do
	printf "\n        Trait: ${trait}"
    	## First, we're going to do the sorting and subsetting common to all
        ## sets so that we don't keep repeating operations that only need to
        ## be done once.
        printf "\n                Sorting and subsetting top ${topn} SNP"
        printf "\n                (Note: There may be more than 1 SNP per gene)"
        for breed in $breeds; do
		# Do the sorting
                sort -g -k 4 "all_snp_${breed}_${trait}_april.txt" -o "all_snp_${breed}_${trait}_april.txt"
                sort -g -r -k 4 "all_snp_${breed}_${trait}_pers_april.txt" -o "all_snp_${breed}_${trait}_pers_april.txt"
                # Keep the top ${topn} from each set
                head -n ${topn} "all_snp_${breed}_${trait}_april.txt" > "all_snp_${breed}_${trait}_april.txt.${topn}"
                head -n ${topn} "all_snp_${breed}_${trait}_pers_april.txt" > "all_snp_${breed}_${trait}_pers_april.txt.${topn}"
                # Sort the subsets by gene ID
                sort -k 1 "all_snp_${breed}_${trait}_april.txt.${topn}" -o "all_snp_${breed}_${trait}_april.txt.${topn}"
                sort -k 1 "all_snp_${breed}_${trait}_pers_april.txt.${topn}" -o "all_snp_${breed}_${trait}_pers_april.txt.${topn}"
                # Compute the intersections
                printf "\n                ${breed} -- ${trait}"
                join -j1 -o 1.2,1.3,1.4,1.5,2.5 <(<"all_snp_${breed}_${trait}_april.txt.${topn}" awk '{print $1"-"$2"-"$3" "$0}' | sort -k1,1) <(<"all_snp_${breed}_${trait}_pers_april.txt.${topn}" awk '{print $1"-"$2"-"$3" "$0}' | sort -k1,1) > "all_snp_${breed}_lh_${trait}_shared.txt"
	        printf "\n                        $(wc -l "all_snp_${breed}_lh_${trait}_shared.txt") genes"
	        # Clean-up
	        rm "all_snp_${breed}_${trait}_april.txt.${topn}"
	        rm "all_snp_${breed}_${trait}_pers_april.txt.${topn}"
        done
done
printf "\n"
