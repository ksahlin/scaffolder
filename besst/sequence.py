class SubSequence(object):
	"""Docstring for Subsequence"""
	reverse_table = {'A':'T','C':'G'}
	def __init__(self, contig_object,contig_start,contig_end):
		super (SubSequence, self).__init__()

		# Contig whereabouts, indexed like lists, i.e.
		# 0-indexed and [:7] means positions 0 to 6
		self.contig = contig_object
		self.contig_start_pos = contig_start
		self.contig_end_pos = contig_end

		# Scaffold whereabouts, indexed like lists, i.e.
		# 0-indexed and [:7] means positions 0 to 6
		self.scaffold_revcomp = None			
		self.scaffold = None
		self.scaffold_start_pos = None
		self.scaffold_end_pos = None

		print 'subseq object created'
		print self.contig_start_pos, self.contig_end_pos
		

	def __str__(self):
		if self.scaffold_revcomp:
			return(reverse_complement(self.contig.sequence[self.contig_start_pos : self.contig_end_pos]))
		else:
			return self.contig.sequence[self.contig_start_pos : self.contig_end_pos]

	def __len__(self):
		return(self.contig_end_pos - self.contig_start_pos)

	def reverse_compelment(self):
		return()

	def __lt__(self,other):
		if self.scaffold_start_pos < other.scaffold_start_pos:
			return(True)
		else:
			return(False)



class Contig(object):
	"""Docstring for Contig"""
	def __init__(self,name,sequence):
		super(Contig, self).__init__()
		self.name = name
		self.sequence = sequence



class Scaffold(object):
	"""docstring for Scaffold"""
	def __init__(self, name):
		super(Scaffold, self).__init__()
		self.name = name
		self.subsequences = set()

	def __str__(self):
		ordered_seqs = self.order_seqs()
		scaf = ''
		current_pos = 0
		for subseq in ordered_seqs:
			gap = subseq.scaffold_start_pos - current_pos
			if gap > 0:
				scaf +='N'*gap
			else:
				scaf +='n'
			scaf += str(subseq)
			current_pos += subseq.scaffold_end_pos

		return(scaf)

	def order_seqs(self):
		return(sorted(self.subsequences))

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
		
