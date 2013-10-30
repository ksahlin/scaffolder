from besst.sequence import Scaffold
from besst import score
from besst import interval
from besst.graph import SequenceGraph as g
from besst.graph import LinearPath
from besst.util import input_
from besst.sequence import SubSequence
from besst.sequence import Contig



##
# Hande input to scaffolding from three different files
# 1. A contig fasta file.
# 2. A contigs to subsequenc file (made when parsing bam file)
# 3. A link file that describes all linkings of contigs (made when parsing bam file). 

contigs ={}
for accession,sequence in input_.get_contigs(open('../test/contig.fa','r')):
	c = Contig(accession,sequence)
	contigs[c.name] = c

subsequences = {}
for subseq_name, contig, start_pos, end_pos in input_.get_subsequences(open('../test/subseqs.txt','r')):
	s =SubSequence( subseq_name, contigs[contig], start_pos, end_pos)
	subsequences[subseq_name] = s
	print subseq_name

G = g()

##
# Initialize nodes
for subseq in subsequences:
	G.add_node(subsequences[ subseq ])


##
# Create edges
for seq1,orientation1, seq2, orientation2, link_count in  input_.get_links(open('../test/links.txt','r')):
	#TODO: Calculate distance here
	G.add_edge((subsequences[ seq1 ] ,orientation1),(subsequences[ seq2 ],orientation2),d=0,s=score.nr_links(link_count))


print len(G.edges())
print G.nodes()

G.remove_self_links()


score_list=[]
for node in G.nodes_iter():
	nbrs = G.neighbors(node)
	if nbrs:
		wip = interval.WeightedIntervalProblem(node)
		for nbr in nbrs:
			i = interval.Interval(start=G[node][nbr]['d'], 
				end=G[node][nbr]['d'] + len(nbr[0]),weight=G[node][nbr]['s'],name=nbr)
			wip.add_interval(i)

		wip.weighted_interval_scheduling()
		score_list.append(wip)

print score_list
print sorted(score_list,reverse=True,key=lambda x:x.score)
visited = set()
for wisp_instance in sorted(score_list,reverse=True,key=lambda x : x.score):
	potential_nodes_to_join =set()
	potential_nodes_to_join.add(wisp_instance.startnode)
	if len(wisp_instance.optimal_path) > 1:
		for seq_obj in map(lambda x: x[3][0], wisp_instance.optimal_path[:-1]):
			potential_nodes_to_join.add(seq_obj,True)
			potential_nodes_to_join.add(seq_obj,False)
	potential_nodes_to_join.add(wisp_instance.optimal_path[-1].name)

	if not visited.intersection(potential_nodes_to_join):
		##
		# make trusted path create an unbreakable linear path in the
		# scaffold graph and returns all the visited nodes
		visited_nodes = G.remove_deactivated_edges(wisp_instance) 
		G.construct_trusted_edges(wisp_instance)

		visited.update(visited_nodes)


print 'GRAPH:'
for edge in G.edges():
	print (edge[0][0]),
	print(edge[1][0])


##
# Make scaffolds

# find start nodes
start_nodes = set()
for node in G.nodes_iter():
	##
	# fount start node
	if not G.neighbors((node[0],True)):
		start_nodes.add((node[0],True))
	if not G.neighbors((node[0],False)):
		start_nodes.add((node[0],False))

print 'start nodes',start_nodes

for node in G.nodes():
	print repr(node[0]), node

visited =set()
scaffold_index = 1
for start_node in start_nodes:
	if start_node not in visited:
		print 'Making scaffold',scaffold_index
		s = Scaffold(scaffold_index)
		path = LinearPath(G,start_node)
		for node,gap in path:
			print node[0].contig.name
			if node[1]:
				node[0].rc = False
				s(str(node[0]))
			else:
				node[0].rc = True
				s(str(node[0]))
			if gap <= 0:
				s('n')
			else:
				s('N'*gap)
		visited.add(start_node)
		visited.add(node)
		print visited
		print '>scaffold'+str(scaffold_index)
		print s
		scaffold_index += 1


