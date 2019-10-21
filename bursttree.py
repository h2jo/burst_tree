import random

def seq_to_dist(seq):
    """Creates an empirical distribution of elements given a sequence.

    Distribution is represented as a dictionary where key is an element and value is the 
    number of such elements.
    """
    dist={}
    for item in seq:
        dist[item]=dist.get(item,0)+1
    return dist

def dist_to_seq(dist):
    """Creates a sequence with arbitrary order of elements from an empirical distribution of elements.
    """
    seq=[]
    for item,num in dist.items():
        seq.extend(num*[item])
    seq.sort(reverse=True)
    return seq

class BurstTree(object):
    def __init__(self):
        """Construct an empty burst tree."""
        #We only store the splits in the tree and no information on the IETs
        #Each split is a tuple of form (left,right) where left is the index of the 
        #split on the left and right is the index of the split on the right. If left
        #or right bracnh is not a split but a leaf then they are set to None.
        #For example:
        # splits = [(1,None),(None,None)] would correspond to  the following tree:
        #         0
        #       /   \
        #      1    leaf
        #    /   \
        # leaf   leaf
        self.splits=[]

    def decompose(self,timeseq):
        """Decomposes the time sequence into burst tree and IET distribution."""
        assert timeseq==sorted(timeseq), "Time sequence needs to be sorted."

        #Find the IETs and build the IET distribution
        iets=[timeseq[i+1]-timeseq[i] for i in range(len(timeseq)-1)]
        iet_dist=seq_to_dist(iets)
        niets=len(iets)

        #Clusters are indexed with their rank-order of the corresponding IETs. Large IET means small index
        #We build a mapping from cluster index to iet index and another one from iet index to cluster index
        cluster_to_iet_index = list(map(lambda x:x[0],sorted(enumerate(iets), key= lambda x:-x[1])))
        iet_index_to_cluster=[None]*niets
        for i in range(niets):
            iet_index_to_cluster[cluster_to_iet_index[i]]=i

        #From now on we do not need the iet values
        del iets

        #For merging clusters we need to know the IETs at their borders
        #We make a dict where key is the cluster index and value is tuple giving the IET indices at the 
        #left and right borders of the cluster
        cluster_to_left_right_iet={}

        #Go through the clusters in reverse order, i.e., the low IET ones first
        for cluster in reversed(range(niets)):
            iet_index = cluster_to_iet_index[cluster]
            
            #Left side of the cluster
            if iet_index>0 and iet_index_to_cluster[iet_index-1]>cluster: #merge
                cluster_left=iet_index_to_cluster[iet_index-1]
                #update the cluster of the left-most iet
                if cluster_left in cluster_to_left_right_iet:
                    iet_left=cluster_to_left_right_iet[cluster_left][0]
                else:
                    iet_left=iet_index-1
                iet_index_to_cluster[iet_left]=cluster
            else: #no merge
                cluster_left=None
                iet_left=iet_index

            #Right side of the cluster
            if iet_index<niets-1 and iet_index_to_cluster[iet_index+1]>cluster: #merge                    
                cluster_right=iet_index_to_cluster[iet_index+1]
                #update the cluster of the right-most iet                                                                                             
                if cluster_right in cluster_to_left_right_iet:
                    iet_right=cluster_to_left_right_iet[cluster_right][1]
                else:
                    iet_right=iet_index+1
                iet_index_to_cluster[iet_right]=cluster
            else: #no merge                                                                                                                  
                cluster_right=None
                iet_right=iet_index

            #update left-right map and the clusters where the IETs belong to
            cluster_to_left_right_iet[cluster]=(iet_left,iet_right)
            iet_index_to_cluster[iet_left]=cluster
            iet_index_to_cluster[iet_right]=cluster

            #Add a new split, we add in reverse order, and reverse the splits at the end
            self.splits.append((cluster_left,cluster_right) )
        
        self.splits.reverse() #reverse the splits

        return iet_dist

    def recompose(self,ietdist,firsttime=0):
        """Recompose this tree structure and the given inter-event time distribution (ietdist).

        Parameters
        ----------
        ietdist : Dict
           The distribution of inter-event times: a dictionary with keys as IETs and values
           as counts of those IETs.
        firsttime: int, optional
           The time of the first event in the recomposed time sequences

        Returns
        -------
        list
           The reconstructed time sequence of the events.
        """
        ietseq=dist_to_seq(ietdist)

        seq=[firsttime]
        for index in self.traverse_from_left():
            iet=ietseq[index]
            seq.append(seq[-1]+iet)
        return seq

    def traverse_from_left(self, next=0):
        """A generator for traversing this tree from the left.
        """
        left = self.splits[next][0]
        right = self.splits[next][1]
        if left!=None:
            yield from self.traverse_from_left(left)
        yield next
        if right != None:
            yield from self.traverse_from_left(right)
            
    def number_of_leaves(self):
        return len(self.splits)+1

    def burst_sizes(self):
        """Returns the burst sizes for all the splits.
        """
        b=[None]*len(self.splits)
        for i in reversed(range(len(self.splits))):
            left=self.splits[i][0]
            right=self.splits[i][1]

            nleft=1 if left==None else b[left]
            nright=1 if right==None else b[right]
            
            b[i]=nleft+nright

        return b


def test_decomposition(size,ietrange=10**8,verbose=True):
    #Generate a random time sequence
    if verbose: print("Generating a random time sequence with", size, "events.")
    seq=sorted(random.sample(range(ietrange),size))

    #Decompose the time sequences
    if verbose: print("Decomposing the time sequence.")
    tree=BurstTree()
    iet_dist=tree.decompose(seq)
    first=seq[0] #save the first event time, as this information is not in the IET or the tree

    #Make sure that we can recompose the time sequence just with the tree, IET distribution and the first IET
    if verbose: print("Recomposing the time sequence.")
    recomposed = tree.recompose(iet_dist,first)
    if verbose: print("Comparing the original time sequence and the recomposed one...")
    assert seq==recomposed
    if verbose: print("Success!")
    

if __name__ == "__main__":
    test_decomposition(10**6)



