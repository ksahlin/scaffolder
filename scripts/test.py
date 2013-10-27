from besst.sequence import SubSequence as ss
from besst.sequence import Scaffold
from besst.sequence import Contig
from besst import score
# we are reading lines from a file with format for each link:
#contig1 	start_cord_cont1 	end_cord_contig1	contig2	 start_cord_cont2  end_cord_cont2 nr_links 
# %observations subseq1 (coordinates on the subseq of contig1 in space separated format with '%' sign first
# %observations subseq1 (coordinates on the subseq of contig1 in space separated format with '%' sign first

# True sequence: AAACCCGGGTTT

##
#preprocessing module

c1 = Contig('contig1','AAAGGG') # this should be split
c2 = Contig('contig2','AAAGGG') # this should be split and rc
print(c1.sequence)
s1 = ss(c1,0,3)
s2 = ss(c1,3,6)
s3 = ss(c2,0,3)
s4 = ss(c2,3,6)
print len(s1),len(s2)

# scaf = Scaffold('scaf1')
# # scaf.add_subsequence(s2,False,0)
# scaf.add_subsequence(s1,False,7)
print(s1)
print(s2)
# print(scaf)




from besst.graph import SequenceGraph as g

G = g()

G.add_node(s1)
G.add_node(s2)
G.add_node(s3)
G.add_node(s4)

G.add_edge((s1,True),(s4,True),d=0,s=score.nr_links(10))
G.add_edge((s2,True),(s3,True),d=0,s=score.nr_links(12))
G.add_edge((s4,False),(s2,False),d=0,s=score.nr_links(7))


#  false link!
G.add_edge((s1,True),(s1,False))

# G.add_edge()
print(G)


for edge in G.iteredges():
	print edge
	print G[edge[0]][edge[1]]['s']

for node in G.


print (G.neighbors((s1,True)))

#print G.most_likely_neighbor((s1,True))[0].contig.name


