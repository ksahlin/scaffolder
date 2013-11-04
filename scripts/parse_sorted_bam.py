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

def main(infile):
	bam_file = bamio.open_bam_file(infile)
	contigs = create_contig_objects(bam_file)

	base_pair_tracker = BasePairTracker(contigs[bam_file.gettid(bam_file.references[0])])
	for read in bam_file:
		if base_pair_tracker.accept(read):
			if read.tid != base_pair_tracker.contig.tid:
				print read.tid, base_pair_tracker.contig.name
				base_pair_tracker.return_stats()
				for bp in base_pair_tracker.contig.breakpoints:
					if bp[0] > 2000 and bp[0] < base_pair_tracker.contig.length - 2000: 
						print base_pair_tracker.contig.breakpoints, 'lenthg',base_pair_tracker.contig.length
				base_pair_tracker = BasePairTracker(contigs[read.tid])

			base_pair_tracker.visit(read)

	base_pair_tracker.return_stats()







if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="Find regions that seem to contain a misassembly.")
    arg_parser.add_argument("alignment_file", type=str, help="The aligned reads to infer contig breakpoints from from.")
    # arg_parser.add_argument("output_file", type=argparse.FileType("w"), help="Output file containing predicted SV regions.")
    args = arg_parser.parse_args()
    main(args.alignment_file)