from collections import defaultdict

##
# This class holds statistics of reads mapped to a reference
# genome. The idea is to feed it aligned reads and for each
# read it is fed, the statistics are updated. You may then at
# any point get the statistics you want.
#
class BasePairTracker(object):
    ##
    # Constructor.
    #
    # @param handle_bp Takes care of a BasePairStats object that have no more 
    #                        stats to update. Takes a tid and BasePairStats object.
    #
    def __init__(self, contig,interval = 100):
        # Stats for each active base pair in their respective chromosome
        self.active_bps =  defaultdict(BasePairStats)
        self.contig = contig
        self.interval = interval

    ##
    # Determines whether we should visit a read or not. In
    # this case we ignore reads that have not been mapped.
    # 
    # @param read The read to visit.
    # 
    # @return True if the read is mapped, false otherwise.
    #
    def accept(self, read):
        #TODO: rewrite this and include consistent orientation detection
        if read.rname != read.mrnm or read.is_unmapped or read.mpos <= read.pos or abs(read.mpos - read.pos) > 6500 or abs(read.mpos - read.pos) < 2000:
            return False
        else:
            return True

    ##
    # Updates the statistics by being fed a new aligned read.
    #
    # @param read An aligned read.
    #
    def visit(self, read):
        #self.update_stats_for_soft_clipped_reads(read)
        self.update_stats_for_span_coverage(read)



    def update_stats_for_span_coverage(self, read):
        start = int(round(read.pos + self.interval, -2))
        end = int(round(read.mpos, -2))

        for position in range(start, end, self.interval):
            #print position
            self.active_bps[ position ].span_coverage += 1

    def breakpoints(self):
        break_coordinates = []
        start = self.interval
        end = int(round(self.contig.length, -2))
        last_break_point = 0
        for position in range(start, end, self.interval)[1:]:
            # print position, next_position
            #TODO: find fancy way of calling a breakpoint here
            # If position does not exist in "active_bps", it returns 0. Thus, it will also enter the if-statement
            if self.active_bps[ position ].span_coverage < 1 :
                #print position - prev_position
                if  position - last_break_point > self.interval:
                    break_coordinates.append( (position, self.active_bps[ position ].span_coverage))

                last_break_point = position
                
        return(break_coordinates)

    # def breakpoints(self):
    #     for point in self.tracked_points:
    #         pass


##
# This is a simple class for holding different sufficient
# statistics of base pairs in the reference genome.
#
class BasePairStats(object):
    def __init__(self):
        # TODO: Eventually include more stats per base pair here
        # # Number of reads ending at the given position, split reads
        # # are considered ending in multiple positions.
        # self.num_ending = 0.0

        # # Number of reads starting at the given position
        # self.num_starting = 0.0

        self.span_coverage = 0


