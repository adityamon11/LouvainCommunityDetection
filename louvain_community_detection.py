"""
louvain_cd.py: Python code for community detection using Louvain Method.

@Author: "Aditya Kulkarni, Ameya Marathe"
"""

import igraph
from collections import defaultdict
from time import time

nodes = []
edges = []
m = 0                                                                           # sum of all of the edge weights in ntwk
node_weights = []                                                               # sum of all weights incedent on each node
node_edges = defaultdict(list)                                                  # dictionary containing all edges incedent on each node
comms = []                                                                      # array of assigned communities to each node
partitions = []                                                                 # actual partitions
sum_in_comm = []
sum_total = []
weights = []

def read_ip(loc):
    """
    Method to process input data passed as .edges file using igraph module.
    Storing it in variables titled nodes, edges.
    nodes - sorted array of nodes
    edges - array of edges in the format [((source, target), weight)]
    Input: location of .txt file
    Output: None
    """
    nodes_tup = []
    global nodes, edges, node_weights
    graph = igraph.read(loc, format="ncol", directed=True)
    print(f"The number of nodes are {len(graph.vs)} and edges are {len(graph.es)}")
    nodes = sorted([i.index for i in graph.vs])
    nodes_tup = [(i.index, i.index) for i in graph.vs]
    nodes_tup.sort(key=lambda x: x[0])
    edges = [((nodes_tup[i.source][1], nodes_tup[i.target][1]), 1) for i in graph.es]

def init_variables():
    """
    Method to initialize and create some base variables.
    Input: None
    Output: None
    """
    global m, node_weights, node_edges, comms, weights
    node_weights = [0] * len(nodes)
    comms = [i for i in nodes]
    weights = [0] * len(nodes)
    for edge in edges:
        m += edge[1]
        node_weights[edge[0][0]] += edge[1]
        node_weights[edge[0][1]] += edge[1]
        node_edges[edge[0][0]].append(edge)
        node_edges[edge[0][1]].append(edge)

def start_louvain():
    """
    Base method to start the execution from.
    Input: None
    Output: Partitions and the modularity score of network
    """
    global partitions
    cur_nodes = [i for i in nodes]
    cur_edges = [i for i in edges]
    best_mod = float("-inf")
    while True:
        partition = modularity_optimization(cur_nodes, cur_edges)
        mod_score = calc_mod(partition)
        partition = [i for i in partition if i]
        if partitions:
            real = []
            for part in partition:
                divisions = []
                for division in part:
                    divisions += partitions[division]
                real.append(divisions)
            partitions = real
        else:
            partitions = partition
        if mod_score == best_mod:
            break
        cur_nodes, cur_edges = community_aggregation(cur_edges, partition)
        best_mod = mod_score
    return (partitions, best_mod)
        
def calc_mod(partition):
    """
    Method to calculate the modularity score for network.
    Input: partitions
    Output: modularity score
    """
    modularity = 0
    m2 = m * 2
    for i in range(len(partition)):
        modularity += sum_in_comm[i] / m2 - (sum_total[i] / m2) ** 2
    return modularity

def calc_mod_gain(node, comm, sum_weight_in_comm):
    """
    Method to calculate the gain for a particular vertex
    Input: vertex, particular community, sum of the weights in community
    Output: modularity gain for a vertex
    """
    return 2 * sum_weight_in_comm - sum_total[comm] * node_weights[node] / m

def modularity_optimization(cur_nodes, cur_edges):
    """
    Method for first phase of Louvain Algorithm, which is modularity optimization.
    Input: lists of nodes and edges.
    Output: current best partitions
    """
    global sum_in_comm, sum_total, comms
    current_best = starting_partition(cur_edges, cur_nodes)
    while True:
        mod_change = 0
        for vertex in cur_nodes:
            curr_comm = comms[vertex]
            loc_best = curr_comm
            best_mod_gain = 0
            current_best[curr_comm].remove(vertex)
            best_common_edges = 0
            for edge in node_edges[vertex]:
                if edge[0][0] == edge[0][1]:
                    continue
                if edge[0][0] == vertex and comms[edge[0][1]] == curr_comm or edge[0][1] == vertex and comms[edge[0][0]] == curr_comm:
                    best_common_edges += edge[1]
            sum_in_comm[curr_comm] -= 2* (best_common_edges + weights[vertex])
            sum_total[curr_comm] -= node_weights[vertex]
            comms[vertex] = -1
            loc_communities = {}
            neighbors = get_neighbors(vertex)
            for neighbor in neighbors:
                loc_community = comms[neighbor]
                if loc_community in loc_communities:
                    continue
                loc_communities[loc_community] = 1
                common_edges = 0
                for edge in node_edges[vertex]:
                    if edge[0][0] == edge[0][1]:
                        continue
                    if edge[0][0] == vertex and comms[edge[0][1]] == loc_community or edge[0][1] == vertex and comms[edge[0][0]] == loc_community:
                        common_edges += edge[1]
                mod_gain = calc_mod_gain(vertex, loc_community, common_edges)
                if mod_gain > best_mod_gain:
                    loc_best = loc_community
                    best_mod_gain = mod_gain
                    best_common_edges = common_edges
            current_best[loc_best].append(vertex)
            comms[vertex] = loc_best
            sum_in_comm[loc_best] += 2* (best_common_edges + weights[vertex])
            sum_total[loc_best] += node_weights[vertex]
            if curr_comm != loc_best:
                mod_change = 1
        if not mod_change:
            break
    return current_best


def get_neighbors(vertex):
    """
    Method to get neighbors for any vertex.
    Input: vertex to find neighbors for
    Output: all edges incedent to that vertex
    """
    for edge in node_edges[vertex]:
        if edge[0][0] == edge[0][1]:
            continue
        if edge[0][0] == vertex:
            yield edge[0][1]
        if edge[0][1] == vertex:
            yield edge[0][0]


def starting_partition(cur_edges, cur_nodes):
    """
    Method to get intial partitions.
    Input: list of nodes and edges
    Output: initial partitions of format [[0], [1], [2], [3]]
    """
    global sum_in_comm, sum_total
    partition = []
    sum_in_comm, sum_total = [], []

    for vertex in cur_nodes:
        partition.append([vertex])
        sum_in_comm.append(0)
        sum_total.append(node_weights[vertex])
    for edge in cur_edges:
        if edge[0][0] == edge[0][1]:
            sum_in_comm[edge[0][0]] += edge[1]
            sum_in_comm[edge[0][1]] += edge[1]                                  # Sum of internal edges are doubled for self loops
    return partition

def community_aggregation(cur_edges, partition):
    """
    Method for second phase of Louvain Algorithm, which is community aggregation.
    Input: list of edges, partitions
    Output: lists of nodes and edges
    """
    global comms, node_weights, node_edges, weights
    loc_nodes = [e for e in range(len(partition))]
    loc_communities = []
    d = {}
    i = 0
    for community in comms:
        if community in d:
            loc_communities.append(d[community])
        else:
            d[community] = i
            loc_communities.append(i)
            i += 1
    comms = loc_communities
    loc_edges = defaultdict(int)
    for edge in cur_edges:
        ci = comms[edge[0][0]]
        cj = comms[edge[0][1]]
        loc_edges[(ci,cj)] += edge[1]
    
    loc_edges = [(s,t) for s,t in loc_edges.items()]
    node_weights = [0] * len(loc_nodes)
    node_edges = defaultdict(list)
    weights = [0] * len(loc_nodes)
    for edge in loc_edges:
        node_weights[edge[0][0]] += edge[1]
        node_weights[edge[0][1]] += edge[1]

        if edge[0][0] == edge[0][1]:
            weights[edge[0][0]] += edge[1]
        node_edges[edge[0][0]].append(edge)
        node_edges[edge[0][1]].append(edge)

    comms = [node for node in loc_nodes]
    return(loc_nodes,loc_edges)



if __name__ == '__main__':
    start = time()
    read_ip("datasets/youtube.txt")
    init_variables()
    partitions, modularity = start_louvain()
    print(partitions)
    print("The number of partitions is - ", len(partitions))
    print("The modularity score of the network is - ", modularity)
    print("The running time is ", time() - start)