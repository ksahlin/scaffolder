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
			has been given an optimal solution (using a common DP- algorithm)
			and removes all edges from intervals that was not included in the optimal solution.
		'''
		# visited = set()
		# remove edges from one end (inner end) from start node
		self.remove_nbr_edges(wisp_instance.startnode)
		# visited.add(wisp_instance.startnode)
		# remove edges from both ends from all inner sequences
		# if more than one neighbor in interval
		if len(wisp_instance.optimal_path) > 1:
			for interval in wisp_instance.optimal_path[:-1]:
				seq_obj = interval.node[0]
				self.remove_nbr_edges((seq_obj,True))
				self.remove_nbr_edges((seq_obj,False))
				# visited.add((seq_obj,True))
				# visited.add((seq_obj,False))

		# remove edges from one end (inner end) from end node
		last_node_in_path = wisp_instance.optimal_path[-1].node
		self.remove_nbr_edges(last_node_in_path)
		# visited.add(last_node_in_path)
		# return(visited)

		#TODO: Double check the sanity of this cycle removal. Create test instance for this
		# Eventual cycle removal:
		# if other end of first node connects to other end of last node in path
		# connect to eachother (cycle).
		other_end_first = (wisp_instance.startnode[0], 1-wisp_instance.startnode[1])
		other_end_last = (last_node_in_path[0],1-last_node_in_path[1])
		if other_end_first == other_end_last:
			#print 'LLLLOOOOOLLL'
			self.remove_nbr_edges(other_end_last)

			
	def construct_trusted_edges(self,wisp_instance):
		''' Takes a weighted interval scheduling problem instance that
			has been given an optimal sulution (using a common DP- algorithm)
			and constructs edges between all intervals that were found in the optimal 
			solution.
		'''
		# make first trusted edge
		self.add_edge(wisp_instance.startnode,wisp_instance.optimal_path[0].node,d=wisp_instance.optimal_path[0].start,s=None)
		# print wisp_instance.score
		# the rest 
		for i1, i2 in zip(wisp_instance.optimal_path, wisp_instance.optimal_path[1:]):
			self.add_edge((i1.node[0],not i1.node[1]),i2.node, d=i2.start - i1.end,s=None)


class LinearPath(SequenceGraph):
	"""Iterator over the remaining linear paths in a 
		SequenceGraph object"""
	def __init__(self, graph, node):
		super(LinearPath, self).__init__()
		self.graph = graph
		self.visited = set()
		self.node = (node[0],not node[1])
		self.exit = False
	def __iter__(self):
		return self

	def next(self):
		if self.exit or self.node in self.visited:
			raise StopIteration

		node = self.node
		self.visited.add(node)
		nbr = self.graph.neighbors(self.node)
		if not nbr:
			self.exit = True
			return(node, 0)

		else:
			self.node = (nbr[0][0],not nbr[0][1])
			gap = self.graph[node][nbr[0]]['d']
			return(node,gap)


