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


	def remove_nbr_edges(self,node):
		for nbr in self.neighbors(node):
			self.remove_edge(node,nbr)

	def remove_deactivated_edges(self,wisp_instance):
		''' Takes a weighted interval scheduling problem instance that
			has been given an optimal sulution (using a common DP- algorithm)
			and removes all edges from intervals that was not included in the optimal solution.
		'''
		visited = set()
		# remove edges from one end (inner end) from start node
		self.remove_nbr_edges(wisp_instance.startnode)
		visited.add(wisp_instance.startnode)
		# remove edges from both ends from all inner sequences
		# if more than one neighbor in interval
		if len(wisp_instance.optimal_path) > 1:
			for interval in wisp_instance.intervals[:-1]:
				seq_obj = interval.name
				self.remove_nbr_edges((seq_obj,True))
				self.remove_nbr_edges((seq_obj,False))
				visited.add((seq_obj,True))
				visited.add((seq_obj,False))

		# remove edges from one end (inner end) from end node
		last_node_in_path = wisp_instance.intervals[-1].name
		self.remove_nbr_edges(last_node_in_path)
		visited.add(last_node_in_path)
		return(visited)
			
	def construct_trusted_edges(self,wisp_instance):
		''' Takes a weighted interval scheduling problem instance that
			has been given an optimal sulution (using a common DP- algorithm)
			and constructs edges between all intervals that were found in the optimal 
			solution.
		'''
		# make first trusted edge
		self.add_edge(wisp_instance.startnode,wisp_instance.intervals[0].name,d=wisp_instance.intervals[0].start,s=None)
		# the rest 
		for i1, i2 in zip(wisp_instance.intervals, wisp_instance.intervals[1:]):
			self.add_edge((i1.name[0],not i1.name[1]),i2.name, d=i2.start - i1.end,s=None)


class LinearPath(SequenceGraph):
	"""Iterator over the remaining linear paths in a 
		SequenceGraph object"""
	def __init__(self, graph, node):
		super(LinearPath, self).__init__()
		self.graph = graph
		self.node = (node[0],not node[1])
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
			return(node,gap)

		






	# def make_scaffolds(self):
	# 	return


		
	# def closest_neighbor(self):
	# 	super(SequenceGraph, self).neighbors(self.closest_neighbor)
	# 	pass


