# This is code provided by Francesco Tiezzt at NCSU on 05/08/14, and
# modified by JBC on 10/09/2014.

rm(list=ls()); gc()

source("http://bioconductor.org/biocLite.R")
library("biomaRt")

#Selecting Biomart Database and Dataset
ensembl<-useMart("ensembl")
#Update to Bos taurus genes (UMD3.1)
ensembl <- useDataset("btaurus_gene_ensembl",mart=ensembl)
#Filters: restriction on the query
#listFilters(ensembl)
#Attributes: values queried
#listAttributes(ensembl)

# set working directory
setwd("/Users/jcole/Documents/AIPL/Genomics/Persistency/networks")

get_gene_names <- function(trait_group, breed, trait_name, num_snp, window_size){
	print(paste("trait_group   :", trait_group))
	print(paste("breed         :", breed))
	print(paste("trait_name    :", trait_name)) 
	print(paste("num_snp       :", num_snp))
        print(paste("window_size   :", window_size))

	# read your snp file
	#snp <- read.table("snp.txt", header=T)
        infilename <- paste("top_200_snp_", breed, "_", trait_group, ".txt", sep="")
        print(paste("Reading data from", infilename))
	ho_pers_all <- read.table(infilename, sep=" ", col.names=c("trait", "snp_rank",
		"snp_name", "snp_chrome", "snp_loc", "effect"))
	print(paste(nrow(ho_pers_all), "records read from", infilename))

	#if ( nrow(ho_pers_all) < 10 ) { print(ho_pers_all) }

	# order it so you have the most important snps on top
	ho_pers <- subset(ho_pers_all, trait==trait_name)
	print(paste(nrow(ho_pers), "records read from", infilename, "for trait", trait_name))
	ho_pers <- ho_pers[order(-ho_pers$effect),]
	print(head(ho_pers))
	
	# set what you want to get, I put it into 'head'
	head <- c("start_position", "end_position", "ensembl_gene_id", "wikigene_name", "wikigene_description", "external_gene_id", "go_id")
	
	# then use it as header of a dataframe where I'll put the genes annotated
	genes <- data.frame(matrix(NA, nrow=0, ncol=length(head))); colnames(genes) <- head
		
	# set filters to use in finding genes, I used to map windows, so I had start and end positions. If you need to annotate SNPs, you can use start=end
	filters <- c("chromosome_name", "start", "end")
		
	# this is to explore around the SNP
	# if you put any value other than 0, you'll explore the window around that marker. This value gives the distance from the marker in base pairs
	# I suggest you to keep a big number, look around, and dig out the closest gene later.
	#dist <- 50000
        dist <- window_size	

	# I use to store otuput into a list, which I can easily recall. The genes around the first SNP will be in the first element, and so on.
	annot <- list()
	
	# the loop will look over your SNP file, picking SNPs from top untill a given row, then you can input in the value 'till'
	#till=200
	till <- num_snp	

	print("Retrieving gene information for each SNP.")
	for(i in 1:till){
	  
		# create a list containing chromosome, start and stop position of the window you want to check
		# If you don't put ay distance, then start=end and you find only genes overlapping the SNP
		position <- list(
			as.numeric(ho_pers$snp_chrome[i]), #chromosome
			as.numeric(ho_pers$snp_loc[i])-dist, # starting position
			as.numeric(ho_pers$snp_loc[i])+dist # ending position
		)
		
		# getBM will go online and dig the genes laying in the window previously defined
		annot[[i]]  <- getBM(attributes=head, filters=filters, values=position, mart=ensembl, verbose=FALSE)
		# you can take a look at what is coming out
		#annot[[i]]
		# remove the info list 
		rm(position)
		# I use to put some sleep between a query and the other, not to 'choke' their system
		Sys.sleep(1)
		if ( i %% 50 == 0 ) {
			print(paste(i,"SNP processed."))
		}
	}

	print("Computing distances from SNP to genes.")
	for(i in 1:till){
		# genes come duplicated because they have a different go_id, we get rid of it and delete duplicatess
		annot[[i]]  <- subset(annot[[i]], !duplicated(ensembl_gene_id), select=(-go_id))
		# this will compute distance between the SNP and each of the genes. 
		# Genes that lay over the SNP will have distance=0, otherwise the distance is with the start or end of the gene, whichever is closer
		d <- matrix(NA, nrow(annot[[i]]),3)
		# Some rows are empty?
		if ( nrow(annot[[i]]) > 0 ) {
			for(j in 1:nrow(annot[[i]])){
				# If there is no gene then annot contains NAs. Francesco's code did not catch that case.
				if ( (!is.na(annot[[i]][j,1])) && (!is.na(annot[[i]][j,1])) ) {
					d[j,1] <- ho_pers$snp_loc[i]-annot[[i]][j,1]
					d[j,2] <- ho_pers$snp_loc[i]-annot[[i]][j,2]
					if(d[j,1]>=0 & d[j,2]<=0){d[j,3] <- 0}else{d[j,3] <- d[j,which(d[j,1:2]==min(d[j,1:2]))]}
				}
			}
  			# the distance will be in variable 'distance'
			annot[[i]]$distance <- d[,3]
		}
	}

	print("Selecting closest gene to each SNP.")
	for(i in 1:till){
		# this will sort all genes from the closest to the farthest, so you can see your gene of interest on top
	        if (!is.null(annot[[i]]$distance)) {
			annot[[i]] <- annot[[i]][order(abs(annot[[i]]$distance)),]
		}
	}

	# I use to put everything on excel spreadsheets so that's easier to read. 
	# On the first sheet, I put the SNPs we queried, then I create a sheet per SNP were I put all genes annotated.
	# The closest gene to every SNP will be on top of each sheet 
	print("Writing SNP and gene information to Excel files.")
	library("xlsx")
	print("\tCreating master XLSX file.")
	write.xlsx(ho_pers[order(ho_pers$snp_chrome,ho_pers$snp_loc),], paste(breed, trait_name, "annotation.xlsx", sep="_"),
	        sheetName="SNP list", row.names=F, col.names=T, append=F)
	print("\tWriting information for each SNP.")
	for(i in 1:till){
        	#print(paste("SNP", i))
		result = tryCatch({
			write.xlsx(annot[[i]], paste(breed, trait_name, "annotation.xlsx", sep="_"),
				sheetName=paste("SNP",i,sep='.'), row.names=F, col.names=T, append=T)
			#print(annot[[i]]$ensembl_gene_id[1])
			if ( i == 1 ) {
				print(paste("Writing results to file ", paste(breed, trait_name, "names.txt", sep="_")))
			}
			write(paste(annot[[i]]$ensembl_gene_id[1],
				ho_pers$snp_name[i],
				ho_pers$effect[i]),
				paste(breed, trait_name, "names.txt", sep="_"),
				append=TRUE)
                        if ( i == 1 ) {
                                print(paste("Writing results to file ", paste(breed, trait_name, "ensembl_names_only.txt", sep="_")))
                        }
			write(paste(annot[[i]]$ensembl_gene_id[1]),
				paste(breed, trait_name, "ensembl_names_only.txt", sep="_"),
				append=TRUE)
		}, warning = function(cond) {
    			#print(paste("There was a warning while processing SNP", i))
			#print(cond)
			this_warning <- TRUE
		}, error = function(cond) {
    			#print(paste("There was an error while processing SNP", i))
			#print(cond)
			this_error <- TRUE
		}, finally={
			#print("Continuing...")
			this_finally <- TRUE
		})
	}
}

# Process persistency traits
num_snp <- 200
window_size <- 20000

trait_group <- "pers"
for ( breed in c("ay", "bs", "ho", "je") ){
	for ( trait_name in c("Milk_pers", "Fat_pers", "Protein_pers", "SCS_pers") ){
		get_gene_names(trait_group, breed, trait_name, num_snp, window_size)
	}
}

# Process yield traits
trait_group <- "yld"
for ( breed in c("ay", "bs", "ho", "je") ){
        for ( trait_name in c("Milk", "Fat", "Protein") ){
                get_gene_names(trait_group, breed, trait_name, num_snp, window_size)
        }       
}

# Remember that SCS is in the NM$ group
trait_group <- "nm"
for ( breed in c("ay", "bs", "ho", "je") ){
        for ( trait_name in c("SCS") ){
                get_gene_names(trait_group, breed, trait_name, num_snp, window_size)
        }       
}
