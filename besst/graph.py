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
			False == Left end of sequence 
			Is sequence links to another sequence from 'True' end, 
			it is forward oriented, otherwise it is reverse complemented"""
		super(SequenceGraph, self).add_nodes_from([(node,True),(node,False)])

	def remove_self_links(self):
		for edge in super(SequenceGraph, self).edges():
			if edge[0][0] == edge[1][0]:
				self.remove_edge(edge[0],edge[1])


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



	def make_trusted_edge(self,node,nbr,distance):
		self.add_edge(node, nbr,s=None,d=distance)

	def make_trusted_path(self,interval_object):
		#TODO: Eventually split this function into two functions.
		visited = set()
		# remove edges from one end (inner end) from start node
		self.remove_nbr_edges(interval_object.startnode)
		visited.add(interval_object.startnode)
		# remove edges from both ends from all inner sequences
		# if more than one neighbor in interval
		if len(interval_object.optimal_path) > 1:
			for interval in interval_object.intervals[:-1]:
				seq_obj = interval[3][0]
				self.remove_nbr_edges((seq_obj,True))
				self.remove_nbr_edges((seq_obj,False))
				visited.add((seq_obj,True))
				visited.add((seq_obj,False))

		# remove edges from one end (inner end) from end node
		last_node_in_path = interval_object.intervals[-1][3]
		self.remove_nbr_edges(last_node_in_path)
		visited.add(last_node_in_path)

		# make first trusted edge
		self.make_trusted_edge(interval_object.startnode,interval_object.intervals[0][3],interval_object.intervals[0][0])
		# the rest 
		for i1, i2 in zip(interval_object.intervals, interval_object.intervals[1:]):
			self.make_trusted_edge((i1[3][0],not i1[3][1]),i2[3], i2[0] - i1[1])

		return(visited)
			
	def remove_nbr_edges(self,node):
		for nbr in self.neighbors(node):
			self.remove_edge(node,nbr)

	# def get_scaffold_path(self,start):
	# 	path = []
	# 	# path.append((start[0],not start[1]))
	# 	for node in self.sequence_generator(start):
	# 		path.append(node)
	# 	return path

	# def path(self,node):
	# 	nbrs = self.neighbors((node[0],not node[1]))
	# 	print nbrs[0]
	# 	if nbrs:
	# 		yield (node[0],not node[1])
	# 		path.append((node[0],not node[1]))
	# 		self.get_scaffold_path(nbrs[0])


class LinearPath(SequenceGraph):
	"""Iterator over the remaining linear paths in a 
		SequenceGraph object"""
	def __init__(self, graph, node):
		super(LinearPath, self).__init__()
		self.graph = graph
		self.node = (node[0],not node[1])
		print self.node[0]
		self.exit = False
	def __iter__(self):
		return self

	def next(self):
		if self.exit:
			raise StopIteration

		node = self.node
		nbr = self.graph.neighbors(self.node)
		if not nbr:
			self.exit = True
			return(node, 0)

		else:
			self.node = (nbr[0][0],not nbr[0][1])
			gap = self.graph[node][nbr[0]]['d']
			print 'Gap:', gap
			return(node,gap)

		






	# def make_scaffolds(self):
	# 	return


		
	# def closest_neighbor(self):
	# 	super(SequenceGraph, self).neighbors(self.closest_neighbor)
	# 	pass


