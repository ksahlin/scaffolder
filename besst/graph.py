from networkx import Graph
import score

class SequenceGraph(Graph):
	"""docstring for Graph"""
	def __init__(self):
		super(SequenceGraph, self).__init__()

	def __str__(self):
		graph_repr = ''
		nodes_visited =set()
		for node,data in self.nodes(data=True):
			if node[0] not in nodes_visited: 
				graph_repr += ''.join(map(lambda x: str(x),[node[0].contig.name,' from ',node[0].contig_start_pos,' to ',node[0].contig_end_pos,' ',data,'.\n']))
			nodes_visited.add(node[0])
		return(graph_repr)

	def add_node(self,node):
		""" True == Right end of sequence,
			False == Left end of sequence """
		super(SequenceGraph, self).add_nodes_from([(node,True),(node,False)])


	def iteredges(self):
		'''
			Iterator over edges that links two different sequences (no self-edges connecting two ends of the same sequence)
		'''
		for edge in super(SequenceGraph, self).edges_iter():
			if edge[0][0] != edge[1][0]:
				yield edge

	def neighbors(self,node):
		'''
			Returns neighbors of a node except itself as a neighbor. (If the contig node is linking to its other end)
		'''
		return filter(lambda x: x[0] != node[0],super(SequenceGraph, self).neighbors(node))


	def make_trusted_edge(self,node,nbr):
		return

	def make_trusted_path(self,edges):
		return
	# def most_likely_neighbor(self,node):
	# 	# TODO: More sophisticated method for most likely neighbor
	# 	return sorted(self.neighbors(node),key=lambda nbr: self[node][nbr]['s'])[0]
			




	# def make_scaffolds(self):
	# 	return


		
	# def closest_neighbor(self):
	# 	super(SequenceGraph, self).neighbors(self.closest_neighbor)
	# 	pass

