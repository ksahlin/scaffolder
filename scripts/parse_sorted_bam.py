import argparse
from besst.util import bamio 
from besst.sequence import Contig
from besst.preprocess.base_pair_tracker import BasePairTracker

def create_contig_objects(bam_file):
	contigs = {}
	for ctg, length in zip(bam_file.references, bam_file.lengths):
		c = Contig( ctg, tid = bam_file.gettid(ctg),length = length)
		contigs[c.tid] = c
	return(contigs)

def create_sequence_file(contig, sequence_file,seq_id,start,stop):
	sequence_file.write(contig.name+'\t'+str(seq_id)+'\t'+ str(start)+'\t'+str(stop)+'\n')
	pass

def main(args):
	bam_file = bamio.open_bam_file(args.alignment_file)
	contigs = create_contig_objects(bam_file)

	base_pair_tracker = BasePairTracker(contigs[bam_file.gettid(bam_file.references[0])])
	sequence_identifier = 1
	for read in bam_file:
		if base_pair_tracker.accept(read):
			if read.tid != base_pair_tracker.contig.tid:
				#print read.tid, base_pair_tracker.contig.name
				if base_pair_tracker.contig.name == 'velvet.93.10':
					print 'HEHEHEHEHEHE'
					print read.tid, base_pair_tracker.contig.name
				breaks = base_pair_tracker.breakpoints()
				# add first and last position of contigs to breaks
				breaks.insert(0,(0,0))
				breaks.append((base_pair_tracker.contig.length,0))
				for start,stop in zip(breaks[:-1],breaks[1:]):
					if base_pair_tracker.contig.length > 4200:  #start[0] > 2000 and bp[0] < base_pair_tracker.contig.length - 2000 and
						if start[0] > 0 or stop[0] < base_pair_tracker.contig.length - 100:
							print start, 'lenthg',base_pair_tracker.contig.length
						create_sequence_file(base_pair_tracker.contig, args.sequence_file,sequence_identifier,start,stop)
						sequence_identifier += 1
				
				base_pair_tracker = BasePairTracker(contigs[read.tid])
			base_pair_tracker.visit(read)
	breaks = base_pair_tracker.breakpoints()







if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="Find regions that seem to contain a misassembly.")
    arg_parser.add_argument("alignment_file", type=str, help="The aligned reads to infer contig breakpoints from from.")
    arg_parser.add_argument("sequence_file", type=argparse.FileType("w"), help="Output sequence file containing\
     all subsequences with positions and what positions witin the contigs.")
    args = arg_parser.parse_args()
    main(args)

