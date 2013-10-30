import os
import pysam

##
# Opens a .bam or .sam file and returns the file
# object.
#
# @param bam_file_path Path to the .bam or .sam.
# 
# @return File object for .bam or .sam.
#
def open_bam_file(bam_file_path):
    bam_file_name, bam_file_ext = os.path.splitext(bam_file_path)
    if bam_file_ext == ".bam":
        return pysam.Samfile(bam_file_path, 'rb')
    elif bam_file_ext == ".sam":
        return pysam.Samfile(bam_file_path, 'r')
    else:
        return IOError("open_bam_file: File must be either .bam or .sam.")

def get_subsequences(subseq_file):
	for line in subseq_file:
		subseq,contig, start_pos, end_pos = line.split()
		yield subseq, contig, int(start_pos), int(end_pos)

def get_links(link_file):
	for line in link_file:
		if line[0] == '%':
			continue
		else:
			s1, orientation1, s2, orientation2, nr_links = line.split()
			o1 = True if int(orientation1) == 1 else False
			o2 = True if int(orientation2) == 1 else False
			yield s1, o1, s2, o2, int(nr_links)


def get_contigs(contig_file):
    k = 0
    sequence = ''
    accession = ''
    for line in contig_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            k += 1
        elif line[0] == '>':
        	yield accession, sequence
        	sequence = ''
        	accession = line[1:].strip().split()[0]
        else:
            sequence += line.strip()
    	yield accession, sequence