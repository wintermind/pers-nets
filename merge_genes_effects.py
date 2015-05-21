#!/Users/jcole/anaconda/bin/python

# We're going to read two files, a BED file and a
# data file, merge the contents on Illumina SNP
# ID, and write an output file for each breed and
# trait.

# Return a dictionary of the traits in each group.
def get_traits():
	traits = {}
	traits['nm'] = ['SCS']
	traits['pers'] = ['Milk_pers', 'Fat_pers', 'Protein_pers', 'SCS_pers']
	traits['yld'] = ['Milk', 'Fat', 'Protein']
	return traits

# Return a dictionary that has Illumina SNP IDs as keys and
# Ensembl gene names as values.
def load_bed(bedfile):
	# The BED file looks like this:
	#
	# ba052234:networks jcole$ head HO_snp_within_25kb_autosomes.txt
	# chr1 120183 BovineHD0100000035 ENSBTAT00000002299
	# chr1 135098 Hapmap43437-BTA-101873 ENSBTAT00000002299
	# chr1 158820 BovineHD0100000048 ENSBTAT00000002299

	#
	# jcole@jcole ~/Documents/AIPL/Genomics/Persistency/networks $ head HO_snp_within_25kb_april.txt
	# chr1 120183 BovineHD0100000035 ENSBTAG00000001753 protein_coding
	# chr1 135098 Hapmap43437-BTA-101873 ENSBTAG00000001753 protein_coding
	# chr1 158820 BovineHD0100000048 ENSBTAG00000001753 protein_coding
	#
	# I'm going to read each line, split it into tokens, and
	# add the gene name (the 4th column) and gene symbol (the
	# 5th column) to a dictionary keyed on SNP name (the 3rd column).
	bed_dict = {}
	inf = open(bedfile, 'r')
	for line in inf:
		pieces = line.split()
		bed_dict[pieces[2]] = "%s %s" % (pieces[3], pieces[4])
	return bed_dict

# Read the files that contain all traits for a group and split
# them into new files, one per trait per breed. The new file is
# space-deliminted and contains three columns: trait name,
# Illumina SNP ID, and the effect size in additive genetic SD.
def split_traits(traitfile, traitgroup, traitname, traitbreed, beddict):

	outfile = 'all_snp_%s_%s_april.txt' % ( traitbreed, traitname )
	outf = open(outfile, 'w')
	inf =  open(traitfile, 'r')
	for line in inf:
		pieces = line.split()
		if pieces[0] == traitname:
			try:
				outline = '%s %s %s\n' % (
					beddict[pieces[2]],
					pieces[2],
					pieces[5] )
				outf.write(outline)
			# The master trait file includes non-
			# autosomal SNP which are not present
			# in the BED file.
			except KeyError:
				pass
	outf.close()

	# The goal is to produce, for each breed and trait, a file
	# that has the following contents:
	#
	# Ensembl gene name  Gene symbol Illumina SNP name     Effect (SD)
	# ENSBTAG00000025200 <foo>       BovineHD1900004545    0.0006783321
	# ENSBTAG00000002728 <bar>       BovineHD0900026933    0.0006514135
	# ENSBTAG00000010818 <arg>       BTA-05125-rs29019289  0.0005757685

if __name__ == '__main__':

	# Get the traits in each group.
	traits = get_traits()

	# Load the BED file.
	bed = load_bed('HO_snp_within_25kb_april.txt')
	#print bed.keys()[0]
	#print bed.values()[0]

	# Create the individual trait files.
	for breed in ['ay', 'bs', 'ho', 'je']:
		for group in ['nm', 'pers', 'yld']:
			traitfile = 'all_snp_%s_%s.txt' % ( breed, group )
			print 'Splitting ', traitfile
			for trait in traits[group]:
				split_traits(traitfile, group, trait, breed, bed)

