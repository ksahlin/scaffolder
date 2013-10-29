from besst.sequence import Scaffold
from besst.sequence import SubSequence as ss
from besst.sequence import Contig
from besst import score
from besst import interval
from besst.graph import SequenceGraph as g
from besst.graph import LinearPath
from besst.util import parse_input

# we are reading lines from a file with format for each link:
#contig1 	start_cord_cont1 	end_cord_contig1	direction	contig2	 start_cord_cont2  end_cord_cont2	direction nr_links 
# %observations subseq1 (coordinates on the subseq of contig1 in space separated format with '%' sign first
# %observations subseq2 (coordinates on the subseq of contig1 in space separated format with '%' sign first

# True sequence: AAACCCGGGTTT

##
#preprocessing module

c1 = Contig('contig1','AAAGGG') # this should be split
c2 = Contig('contig2','AAAGGG') # this should be split and rc

print(c1.sequence)
s1 = ss('s1',c1,0,3)
s2 = ss('s2',c1,3,6)
s3 = ss('s3',c2,0,3)
s4 = ss('s4',c2,3,6)
print len(s1),len(s2)


# scaf = Scaffold('scaf1')
# # scaf.add_subsequence(s2,False,0)
# scaf.add_subsequence(s1,False,7)
print(s1)
print(s2)
# print(scaf)


G = g()

G.add_node(s1)
G.add_node(s2)
G.add_node(s3)
G.add_node(s4)

G.add_edge((s1,True),(s4,True),d=0,s=score.nr_links(10))
G.add_edge((s2,True),(s3,True),d=0,s=score.nr_links(12))
G.add_edge((s4,False),(s2,False),d=0,s=score.nr_links(7))


#  false link!
G.add_edge((s1,True),(s1,False),d=0,s=score.nr_links(5))

# G.add_edge()

print 'GRAPH:'
for edge in G.edges():
	repr(edge[0][0])
	repr(edge[1][0])

print len(G.edges())

G.remove_self_links()
# for edge in G.iteredges():
	# print edge
	# print G[edge[0]][edge[1]]['s']

score_list=[]
for node in G.nodes_iter():
	nbrs = G.neighbors(node)
	if nbrs:
		i = interval.Interval(node)
		for nbr in nbrs:
			start = G[node][nbr]['d']  
			end = G[node][nbr]['d'] + len(nbr[0])
			weight = G[node][nbr]['s']
			i.add_interval(nbr,start,end,weight)
		# print i.intervals

		i.weighted_interval_scheduling()
		score_list.append(i)

print score_list
print sorted(score_list,reverse=True,key=lambda x:x.score)
visited = set()
for interval_object in sorted(score_list,reverse=True,key=lambda x : x.score):
	potential_nodes_to_join =set()
	potential_nodes_to_join.add(interval_object.startnode)
	if len(interval_object.optimal_path) > 1:
		for seq_obj in map(lambda x: x[3][0], interval_object.optimal_path[:-1]):
			potential_nodes_to_join.add(seq_obj,True)
			potential_nodes_to_join.add(seq_obj,False)
	potential_nodes_to_join.add(interval_object.optimal_path[-1][3])

	if not visited.intersection(potential_nodes_to_join):
		# print 'enter here'
		# print len(G.edges())

		# make trusted path create an unbreakable linear path in the
		# scaffold graph and returns all the visited nodes
		visited_nodes = G.make_trusted_path(interval_object) 

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
		scaffold_index += 1
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


print s


# print (G.neighbors((s1,True)))

#print G.most_likely_neighbor((s1,True))[0].contig.name


