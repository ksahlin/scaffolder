import argparse
from besst.util import bamio 
from besst.sequence import Contig
from besst.preprocess.base_pair_tracker import BasePairTracker
from besst.util import input_


def create_contig_objects(bam_file):
	contigs = {}
	for ctg, length in zip(bam_file.references, bam_file.lengths):
		c = Contig( ctg, tid = bam_file.gettid(ctg),length = length)
		contigs[c.tid] = c
	return(contigs)

def create_sequence_file(contig, sequence_file,seq_id,start,stop):
	sequence_file.write(str(seq_id)+'\t'+contig.name+'\t'+ str(start)+'\t'+str(stop)+'\n')

def get_subseqs_of_contigs(args,contig,breaks,sequence_identifier):
	# add first and last position of contigs to breaks
	breaks.insert(0,(0,0))
	breaks.append((contig.length,0))
	for tuple1,tuple2 in zip(breaks[:-1],breaks[1:]):
		start = tuple1[0]
		stop  = tuple2[0]
		if stop - start > 500:
			create_sequence_file(contig, args.sequence_file, sequence_identifier,start,stop)
			contig.add_subsequence(start,stop,sequence_identifier)
			sequence_identifier += 1


	return(sequence_identifier)

def main(args):
	bam_file = bamio.open_bam_file(args.alignment_file)
	contigs = create_contig_objects(bam_file)

	base_pair_tracker = BasePairTracker(contigs[bam_file.gettid(bam_file.references[0])])
	sequence_identifier = 1
	visited_contigs = set()
	for read in bam_file:
		if read.is_unmapped:
			continue
		##
		# For the preprocessing module. Here we split up contigs based
		# on too little support on spanning matepairs.
		if contigs[read.tid].length > 4200:
			if base_pair_tracker.accept(read):
				if read.tid != base_pair_tracker.contig.tid:
					breaks = base_pair_tracker.breakpoints()
					sequence_identifier = get_subseqs_of_contigs(args, base_pair_tracker.contig, breaks, sequence_identifier)			
					base_pair_tracker = BasePairTracker(contigs[read.tid])
				base_pair_tracker.visit(read)
		elif read.tid not in visited_contigs:
			visited_contigs.add(read.tid)
			create_sequence_file(contigs[read.tid], args.sequence_file, sequence_identifier,0,contigs[read.tid].length)
			contigs[read.tid].add_subsequence(0,contigs[read.tid].length,sequence_identifier)
			sequence_identifier += 1

	# last contig to split into subsequences
	breaks = base_pair_tracker.breakpoints()
	sequence_identifier = get_subseqs_of_contigs(args, base_pair_tracker.contig, breaks, sequence_identifier)
	bam_file.close()
	args.sequence_file.close()


	#TODO: We need to check if read is aligned forward or reverse to contig to get correct orientation!
	bam_file = bamio.open_bam_file(args.alignment_file)
	curr_seq = 0
	prev_contig = 0
	connection_handler = input_.SequenceConnectionHandler(curr_seq)
	visited_contigs = set() 
	for read in bam_file:
		if read.tid != read.rnext and not read.is_unmapped and not read.mate_is_unmapped and not read.rnext in visited_contigs:
			visited_contigs.add(read.tid)
			pos1,seq1,len1 = contigs[read.tid].get_subseq_pos_len(read.pos)
			pos2,seq2,len2 = contigs[read.rnext].get_subseq_pos_len(read.pnext)
			if seq1 != curr_seq:
				if curr_seq != 0:
					args.link_file.write(str(connection_handler))
				connection_handler = input_.SequenceConnectionHandler(seq1)
				curr_seq = seq1
			#print pos2,seq2,contigs[read.tid].subsequences,contigs[read.rnext].subsequences
			if seq1 and seq2:
				connection_handler.add_link(len1, pos1, read.is_reverse, seq2, len2, pos2, read.mate_is_reverse)
				#print seq1, pos1, read.is_reverse, seq2 , pos2, read.mate_is_reverse
				#args.link_file.write(str(seq1)+'\t'+contig.name+'\t'+ str(start)+'\t'+str(stop)+'\n')
	#print str(connection_handler)







if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="Find regions that seem to contain a misassembly.")
    arg_parser.add_argument("alignment_file", type=str, help="The aligned reads to infer contig breakpoints from from.")
    arg_parser.add_argument("sequence_file", type=argparse.FileType("w"), help="Output sequence file containing\
     all subsequences with positions and what positions witin the contigs.")
    arg_parser.add_argument("link_file", type=argparse.FileType("w"), help="Output sequence link file containing\
     links between subsequences.")
    args = arg_parser.parse_args()
    main(args)

