import argparse

from besst.sequence import Scaffold
from besst import score
from besst import interval
from besst.graph import SequenceGraph as g
from besst.graph import LinearPath
from besst.util import input_
from besst.sequence import SubSequence
from besst.sequence import Contig


def read_input(args):
	'''
	Handle input to scaffolding from three different files
	1. A contig fasta file.
	2. A contigs to subsequenc file (made when parsing bam file)
	3. A link file that describes all linkings of contigs (made when parsing bam file). 
	'''
	contigs ={}
	for accession,sequence in input_.get_contigs(open(args.contigs,'r')):
		c = Contig(accession, sequence = sequence)
		contigs[c.name] = c

	subsequences = {}
	for subseq_name, contig, start_pos, end_pos in input_.get_subsequences(open(args.seqs,'r')):
		s =SubSequence( subseq_name, contigs[contig], start_pos, end_pos)
		subsequences[subseq_name] = s

	#print len(subsequences)
	G = g()

	##
	# Initialize nodes
	for subseq in subsequences:
		G.add_node(subsequences[ subseq ])


	##
	# Create edges

	for (seq1,orientation1, seq2, orientation2, link_count, gap), naive_gap in  input_.get_links(open(args.links,'r')):
		#TODO: Calculate distance here
		#TODO: Add threshold parameter
		#if link_count >= 5:

		G.add_edge((subsequences[ seq1 ] ,orientation1),(subsequences[ seq2 ],orientation2),d=naive_gap,s=score.nr_links(link_count))

	G.remove_self_links()

	return(G)



def compute_weighted_interval_solutions(G,overlap):
	score_list=[]
	for node in G.nodes_iter():
		nbrs = G.neighbors(node)
		if nbrs:
			wip = interval.WeightedIntervalProblem(node)
			for nbr in nbrs:
				i = interval.Interval(start=G[node][nbr]['d'], 
					end=G[node][nbr]['d'] + len(nbr[0]),weight=G[node][nbr]['s'],node=nbr)
				wip.add_interval(i)

			wip.weighted_interval_scheduling(overlap)
			wip.intervals = []
			score_list.append(wip) 

	return(score_list)


def make_trusted_paths(G,score_list):
	visited = set()
	for wisp_instance in sorted(score_list,reverse=True,key=lambda x : x.score):
		if wisp_instance.score > 0: # 0 score if no solution was found, e.g. too large overlap
			potential_nodes_to_join = set()
			potential_nodes_to_join.add(wisp_instance.startnode)
			#print wisp_instance.score,wisp_instance.optimal_path[0].node[0].contig.name
			if len(wisp_instance.optimal_path) > 1:
				for seq_obj in map(lambda x: x.node[0], wisp_instance.optimal_path[:-1]):
					potential_nodes_to_join.add((seq_obj,True))
					potential_nodes_to_join.add((seq_obj,False))
			potential_nodes_to_join.add(wisp_instance.optimal_path[-1].node)

			if not visited.intersection(potential_nodes_to_join):
				##
				# make trusted path create an unbreakable linear path in the
				# scaffold graph and returns all the visited nodes
				G.remove_deactivated_edges(wisp_instance) 
				G.construct_trusted_edges(wisp_instance)

				visited.update(potential_nodes_to_join)
		else:
			G.remove_deactivated_edges(wisp_instance) 


def make_scaffolds(G,scaffold_index):
	##
	# Make scaffolds

	# TODO: Take care of cycles
	# Two different types: those which a start node and those without start node

	# find start nodes
	start_nodes = set()
	for node in G.nodes_iter():
		##
		# fount start node
		if not G.neighbors(node):
			start_nodes.add(node)


	#print len(G.nodes())
	#print len(start_nodes)
	#TODO: make visited some container in class SequenceGraph
	visited =set()
	
	for start_node in start_nodes:
		#print start_node
		if start_node not in visited:
			s = Scaffold(scaffold_index)
			path = LinearPath(G,start_node)
			for node,gap in path:
				#print node, node[0].name,node[0].contig.name #, node[0].contig.sequence
				if node[1]:
					node[0].rc = False
					s(str(node[0]))
				else:
					node[0].rc = True
					s(str(node[0]))
				if gap <= 0:
					#TODO: Eventually merge/remove overlapping sequences
					s('n')
				else:
					s('N'*gap)

				#add all nodes that has been visited
				visited.add(node)
				visited.add((node[0],1-node[1]))

			visited.add(start_node)

			print '>scaffold'+str(scaffold_index)
			print s
			scaffold_index += 1

	G.remove_nodes_from(visited)
	return(scaffold_index)

def remove_cycles(G):
	# Remove visited nodes from G and check if there are
	# any cyclic structures left
	#for node in G.nodes():
	#	print node[0].contig.name, 'nr nbrs: ', G.neighbors(node), len(G.neighbors(node))
	#print 'cyc',len(G.nodes())

	# All cyclic structures left remove one edge
	# in each cycle in a smart way
	# TODO: Temporary solution!!
	rdm_edge = G.edges()[0]
	#print rdm_edge
	G.remove_edge(rdm_edge[0],rdm_edge[1])

def remove_non_linear_edges(G):
	for node in G.nodes():
		if len(G.neighbors(node)) > 1:
			G.remove_nbr_edges(node)
			#print 'removed edges from ', node
			#print node in G
			#print (node[0],1- node[1]),len(G.neighbors((node[0],1- node[1])))

def main(args):
	G = read_input(args)
	score_list = compute_weighted_interval_solutions(G,args.overlap)
	make_trusted_paths(G,score_list)
	# We now expect only linear scaffolds in a perfect world.
	# We have to remove nodes left with multiple neighbours
	# TODO: create test instance for this!
	remove_non_linear_edges(G)
	scaffold_index = 1
	scaffold_index = make_scaffolds(G,scaffold_index)

	#remove all cycles by removing one edge at a time
	while G.nodes():
		#print 'loop',len(G.nodes())
		remove_cycles(G)
		scaffold_index = make_scaffolds(G,scaffold_index)


if __name__ == "__main__":

    ##
    # TODO: Temporary input reader. Change when preprocedd module is implemented.
    parser = argparse.ArgumentParser(description="besst2")
    parser.add_argument('contigs', type=str, help='contigs fasta file')
    parser.add_argument('links', type=str, help='Links txt file.')
    parser.add_argument('seqs', type=str, help='Sequence file.')
    parser.add_argument('overlap', type=int, help='Overlap parameter')
    args = parser.parse_args()
    main(args)


