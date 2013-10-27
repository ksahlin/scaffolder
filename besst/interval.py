import collections 
import bisect

"""
This algorithm returns the weight of the maximally weighted set of independent intervals
in an interval graph.  Commonly used for scheduling optimization problems.

The algorithm works as follows:

    Imagine this scenario:  You;ve travelling to SV to raise money for you fledging startup,
    and need to schedule investor meetings during your stay to maximize the money you will potentially raise.
    Every potential investor meeting you can take has a start and end time, and you can only be at
    one meeting at a time.  In addition to start and end times, each meeting has value associated with it
    that corresponds to the expected value of money you will raise by taking the meeting.  

    Here's a trivial case (notice that the weight of the interval is independent of its length.):


      *----- 500K -----*    *--- 600K ---*       *-------- 100K --------*
    12pm              1pm  2pm         2:30pm   4pm                    8pm


    In this trival case, the maximall weighting is simply the sum of all weights, since all intervals are independent
    and can be scheduled together.

    In non-trivial cases, the intervals will overlap and tough decisions need to be made.
      
The algorithm's runtime is O(n log n) in the worst case
"""

class Interval(object):
    """docstring for Interval"""
    def __init__(self,startnode):
        super(Interval, self).__init__()
        self.startnode = startnode
        self.score = None
        self.optimal_path = []
        self.intervals = []

    def add_interval(self,seq_obj,start,end,weight):
        self.intervals.append((start,end,weight,seq_obj)) 

    def weighted_interval_scheduling(self):
        '''
        Input a graph G, whose structure is a list of tuples, where each tuple represents
        a weight interval with structure (<start>, <end>, <weight>), where all components are
        less-than comparable types

        --- Doctest ---

        >>> G = [(43,70,27),(3,18,24),(65,99,45),(20,39,26),(45,74,26),(10,28,20),(78,97,23),(0,9,22)]
        >>> i = Interval('startnode')
        >>> for k,x in enumerate(G): i.add_interval('obj'+str(k),x[0],x[1],x[2])
        >>> i.weighted_interval_scheduling()
        100
        >>> G = [(1,5,10),(6,10,12),(1,10,15)]
        >>> i = Interval('startnode')
        >>> for k,x in enumerate(G): i.add_interval('obj'+str(k),x[0],x[1],x[2])
        >>> i.weighted_interval_scheduling()
        22
        '''
        G = self.intervals
        S = collections.defaultdict(int)
        G.sort(lambda x,y: x[1]-y[1])

        start = [x[0] for x in G]  
        end = [x[1] for x in G]

        S[0] = 0
        S[1] = G[0][2]
        for i in xrange(2, len(G)+1):
            S[i] = max(S[i-1], G[i-1][2] + S[ bisect.bisect_left(end, start[i-1]) ])

        self.score = S[len(G)]

        def find_optimal_path(j):
            """ Recursive function to track down the path that gave rise to the optimal solution
            """
            if j == 0:
                return
            else:
                if G[j-1][2] + S[ bisect.bisect_left(end, start[j-1]) ] >= S[j-1]:
                    find_optimal_path(bisect.bisect_left(end, start[j-1]))
                    self.optimal_path.append(self.intervals[j-1])
                else:
                    return find_optimal_path(j-1)

        find_optimal_path(len(G))
        return S[len(G)]



if __name__ == "__main__":
    import doctest
    doctest.testmod()