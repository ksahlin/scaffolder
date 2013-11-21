import os
import sys

import pysam

from collections import defaultdict
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

def get_orientation(o,s1,s2):
    if int(o) == 1:
        return True 
    elif int(o) == 0:
        return False
    else:
        sys.stderr.write('Orientation incorrectly specified for link: {0}, {1}.'.format(s1,s2))
        raise IOError

def get_links(link_file):
	for line in link_file:
		if line[0] == '#':
			continue
		else:
			s1, o1, s2, o2, nr_links, gap = line.split()
			yield s1, get_orientation(o1,s1,s2), s2, get_orientation(o2,s1,s2), int(nr_links), int(gap)


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

class SequenceConnectionHandler(object):
    """docstring for SequenceConnectionHandler"""
    def __init__(self, id_nr):
        super(SequenceConnectionHandler, self).__init__()
        self.id_nr = id_nr
        self.connections_forward = defaultdict(list)
        self.connections_reverse = defaultdict(list)        
    def add_link(self,len1,pos1,is_rc1,other,len2,pos2,is_rc2):
        #TODO: change 100 to parameter.readlength
        if is_rc1:
            obs1 = pos1 + 100
        else:
            obs1 = len1 - pos1
        if is_rc2:
            obs2 = pos2 + 100
        else:
            obs2 = len2 - pos2

        if not is_rc1:
            self.connections_forward[(other,1-is_rc2)].append((obs1,obs2))
        else:
            self.connections_reverse[(other,1-is_rc2)].append((obs1,obs2))

    def make_connections_to_string(self,dict_,not_rc):
        i_am = ''
        for edge in dict_:
            i_am += '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(self.id_nr, not_rc, edge[0], str(edge[1]), len(dict_[edge]),0) #replace 0 with gap
            obs_first = '# '
            obs_second = '# '
            for link in dict_[edge]:
                obs_first += str(link[0]) + ' '
                obs_second += str(link[1]) + ' '
                
            i_am += obs_first+'\n' + obs_second +'\n'
        return(i_am)

    def __str__(self):
        return self.make_connections_to_string(self.connections_forward,1) + self.make_connections_to_string(self.connections_reverse,0)

