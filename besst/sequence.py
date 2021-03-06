import bisect

class SubSequence(object):
	"""Docstring for Subsequence"""
	#TODO: Look up what the reverse complement is for all other charachters Y,R,S etc...
	reverse_table = {'A':'T','C':'G','G':'C','T':'A', 'N':'N',
	'Y':'Y','R':'R','S':'S','M':'M','K':'K','W':'W','n':'n','g':'c','c':'g','t':'a','a':'t'}
	def __init__(self, name, contig_object,contig_start,contig_end):
		super (SubSequence, self).__init__()
		self.name = name
		# Contig whereabouts, indexed like lists, i.e.
		# 0-indexed and [:7] means positions 0 to 6
		self.contig = contig_object
		self.contig_start_pos = contig_start
		self.contig_end_pos = contig_end

		# Scaffold whereabouts, indexed like lists, i.e.
		# 0-indexed and [:7] means positions 0 to 6
		self.rc = None
		self.scaffold = None
		self.scaffold_start_pos = None
		self.scaffold_end_pos = None
		

	def __str__(self):
		if self.rc:
			return self.reverse_compelment()
		else:
			return self.contig.sequence[self.contig_start_pos : self.contig_end_pos]

	def __repr__(self):
		return "From {0}, with coord {1} to {2}".format(self.contig.name,self.contig_start_pos,self.contig_end_pos)

	def __len__(self):
		return(self.contig_end_pos - self.contig_start_pos)

	def __lt__(self,other):
		if self.scaffold_start_pos < other.scaffold_start_pos:
			return(True)
		else:
			return(False)

	def reverse_compelment(self):
		string = self.contig.sequence[self.contig_start_pos : self.contig_end_pos]
		rc_string = ''.join([SubSequence.reverse_table[nucl] for nucl in reversed(string)])
		return(rc_string)




class Contig(object):
	"""Docstring for Contig"""
	def __init__(self,name, tid=None, sequence = None, length = None):
		super(Contig, self).__init__()
		self.name = name
		self.tid = tid
		self.sequence = sequence
		self.length = length
		self.breakpoints = []
		self.neighbors = {}
		self.subsequences = []

	def add_subsequence(self,start,stop,seq):
		self.subsequences.append((start,stop,seq))

	def get_subseq_pos_len(self,contig_pos):
		for start,stop,seq in self.subsequences:
			if contig_pos >= start and contig_pos <= stop:
				return contig_pos - start,seq, stop - start + 1 

		
		return -1,False,-1


class Scaffold(object):
	"""docstring for Scaffold"""
	def __init__(self, name):
		super(Scaffold, self).__init__()
		self.name = name
		self.subsequences = []

	def __call__(self,sequence_string):
		self.subsequences.append(sequence_string)

	def __str__(self):

		return(''.join(self.subsequences))


	def add_subsequence(self,subseq,revcomp,start):
		self.subsequences.add(subseq)
		subseq.scaffold_revcomp = revcomp			
		subseq.scaffold = self.name
		subseq.scaffold_start_pos = start
		subseq.scaffold_end_pos = start + len(subseq)
								

class SplitSequences(object):
	"""docstring for SplitSequences"""
	def __init__(self, arg):
		super(SplitSequences, self).__init__()
		self.arg = arg
		
