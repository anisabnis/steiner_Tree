
''' The prize collecting steiner tree problem.
- Casts the problem as an Integer Program.
- Relaxed to a Linear Program.
- All vertices which have an associated y_i greater than some fraction are considered as terminals.
- A steiner tree is constructed using the primal dual algorithm as mentioned in exercise 7.6.  
'''

import sys
import networkx as nx
from collections import defaultdict
import numpy as np
import itertools
import random
import math
import copy
from scipy.optimize import linprog
from networkx.algorithms.connectivity import local_edge_connectivity
import matplotlib.pyplot as plt

dir = str(sys.argv[1])
gamma = float(sys.argv[2])

def findsubsets(S,m):
    return set(itertools.combinations(S, m))

def seperation_oracle(x, G, root):
    G1 = nx.DiGraph()

    all_edges = list(G.edges)

    for i in range(len(all_edges)):
        edge = all_edges[i]

        G1.add_edge(edge[0], edge[1], capacity=x[i])
        
    all_vertices = list(G.nodes())

    for i,r in enumerate(all_vertices):
        if r == root:
            continue

        flow_value = nx.maximum_flow_value(G1,root,r)


        if flow_value < x[len(all_edges) + i]:
            return [False, r]

    return [True, r]


def parse_input(path):
    ''' Get the penalties for each vertex '''

    ''' Get the links '''
    
    penalties = defaultdict(float)
    edge_weights = defaultdict(float)
    
    G = nx.DiGraph()
   
    f = open(path + '/edges.txt', 'r')
    
    for l in f:
        l = l.strip().split(' ')
        edge = ((int(l[0]), int(l[1])))
        edge_weights[edge] = float(l[2])
        G.add_edge(edge[0],edge[1])
    
    f.close()
        
    f = open(path + '/penalties.txt', 'r')

    for l in f:
        l = l.strip().split(' ')
        penalties[int(l[0])] = float(l[1])
        root = int(l[0])


    f.close()
     
    return [G, edge_weights, penalties, root]



def draw(G):
    pass

# Objective function for the Steiner Tree Problem
def objective(G, edge_weights, penalties):
    c = []

    for e in G.edges:
        c.append(edge_weights[e])

    for v in G.nodes:
        c.append(-1 * penalties[v])
    
    return c            

def constraint_subset(A, B):

    sample_sz = math.ceil(math.log2(len(A)))
    sample = list(range(len(A)))

    sample = random.sample(sample,sample_sz)    

    req_A = np.array([A[i] for i in sample])
    req_B = np.array([B[i] for i in sample])

    return [req_A, req_B, sample]

## Constraints for Steiner Tree Problem
 
def constraints(G, edge_weights, penalties, root):

    all_vertices = copy.deepcopy(list(G.nodes))

    n_vertices = copy.deepcopy(list(G.nodes))
    n_vertices.remove(root)
    all_edges = copy.deepcopy(list(G.edges))
    constraints_by_vertex = defaultdict(list)

    A = np.array([])
    B = []

    c_no = 0

    # Constraints for there exists at least on edge crossing a subset of vertices
    for k in range(1,len(all_vertices)):

        if all_vertices[k] == root:
            continue

        ss = findsubsets(n_vertices, k)

        for s in ss:

            req_edges = [e for e in all_edges if len(set(e).intersection(s)) == 1]

            for i in all_vertices:

                a = [0] * (len(all_edges) + len(all_vertices))
                
                for e in req_edges:
                    a[all_edges.index(e)] = -1

                a[len(all_edges) + all_vertices.index(i)] = 1

                if len(A) == 0:
                    A = np.array(a)
                else:
                    A = np.vstack((A,a))                    

                B.append(0)
                
                constraints_by_vertex[i].append(c_no)

                c_no += 1 
        
    B = np.array(B)

    # All variables greater than zero
    A1 = np.identity(len(all_edges) + len(all_vertices))
    A1 = -1 * A1
    B1 = np.array([0] * (len(all_edges) + len(all_vertices)))

    # All variables less than 1
    A2 = np.identity(len(all_edges) + len(all_vertices))
    B2 = np.array([1] * (len(all_edges) + len(all_vertices)))
    A1 = np.vstack((A1, A2))
    B1 = np.hstack((B1, B2))

    # v_r = 1. The rooted vertex should be 1
    a = [0] * (len(all_edges) + len(all_vertices))
    a[len(all_edges) + all_vertices.index(root)] = 1
    A2 = a
    B2 = np.array([1])

    return [A, A1, A2, B, B1, B2, constraints_by_vertex]


## Returns true if all pairs of end nodes in the st_list are connected in H
def connected_all(st_lst, H):
    for edge in st_lst:
        u = edge[0]
        v = edge[1]

        if local_edge_connectivity(H, u, v)>0:
            continue
        else:
            return False
			
    return True


# Get set of edges which crosses the subsets
def get_subset(C_lst, st_lst):
    C_subset=[]

    for C in C_lst:

        for edge in st_lst:
            s = edge[0]
            t = edge[1]

            if s in C and not t in C and not C in C_subset:
                C_subset.append(C)
            elif t in C and not s in C and not C in C_subset:
                C_subset.append(C)

    return C_subset

# Solve the steiner problem for the selected source-destination pairs
def steiner_primal_dual(G, st_lst, weights):

    H = nx.Graph()

    for u in G.nodes():
        H.add_node(u)

    edgeconstraints = weights
    edge_order = []

    while not connected_all(st_lst, H):

        # Get all connected components
        C_complete = nx.connected_components(H)        
        
        # Get a subset of connected components in which for every s-t i st_list
        # there exists either s or t in the component and not both
        C_subset = get_subset(C_complete, st_lst)
	
        edge_lst=[]	
        
        # Find all boundary edges
        for C in C_subset:
            deltaC = nx.edge_boundary(G,C)
            edge_lst.extend(deltaC)

        edge_lst = list(set(edge_lst))

        min = 10000

        # Find edge with minimum cost to make it tight
        for e in edge_lst:
            u = e[0]
            v = e[1]

            if e in edgeconstraints:                
                if edgeconstraints[((u,v))] < min:
                    min=edgeconstraints[(u,v)]

            elif ((v,u)) in edgeconstraints:
                
                if edgeconstraints[(v,u)] < min:
                    min=edgeconstraints[(v,u)]

        modified = {}

        ## Add the minimum cost edge to the solution
        ## subtract the minimum cost from other edges as well
        ## If the cost for the other edges become tight, add them to the graph as well

        for C in C_subset:

            for edge in nx.edge_boundary(G,C):

                if edge in modified and modified[edge] == 1:
                    H.add_edge(edge[0],edge[1])
                    edge_order.append((edge[0],edge[1]))
                    continue
                else:
                    if edge in edgeconstraints:
                        edgeconstraints[edge] = edgeconstraints[edge] - min

                        if edgeconstraints[edge] <= 0:

                            modified[edge] = 1
                            modified[((edge[1],edge[0]))] = 1
                            H.add_edge(edge[0],edge[1])
                            edge_order.append(((edge[0],edge[1])))

                        else:
                            modified[edge]=0
                            modified[((edge[1],edge[0]))]=0
                    else:
                        edgeconstraints[((edge[1],edge[0]))] = edgeconstraints[(edge[1],edge[0])] - min

                        if edgeconstraints[((edge[1], edge[0]))] <= 0:
                            modified[edge] = 1
                            modified[((edge[1],edge[0]))] = 1
                            H.add_edge(edge[0], edge[1])
                            edge_order.append(((edge[0],edge[1])))
                        else:
                            modified[edge]=0
                            modified[((edge[1],edge[0]))]=0

    iter = len(edge_order) - 1
    
    ## Remove unnecessay edges
    while iter >= 0:

        ed = edge_order[iter]

        if ed in H.edges():
            H.remove_edge(ed[0], ed[1])

        if connected_all(st_lst, H):
            iter -= 1
            continue

        H.add_edge(ed[0], ed[1])
        iter-=1
            
    return H
            

G, edge_weights, penalties, root = parse_input(dir + '/')
[A, A1, A2, B, B1, B2, constraints_by_vertex] = constraints(G, edge_weights, penalties, root)

req_A = A
req_B = B

A2 = np.vstack((A2, A2))
B2 = np.vstack((B2, B2))

A_subset = np.vstack((req_A, A1))
B_subset = np.hstack((req_B, B1))

obj = objective(G, edge_weights, penalties)
sol = linprog(obj, A_ub=A_subset, b_ub=B_subset, A_eq=A2, b_eq=B2, options={"disp": True},method='simplex')

all_vertices = list(G.nodes())
all_edges = list(G.edges())

solution = list(sol.x)

random_nu = 0
while random_nu < gamma:
    random_nu = random.random()
    print(random_nu)

req_vertices = []
for i in range(len(all_vertices)):
    if float(sol.x[len(all_edges) + i]) >= random_nu:
        req_vertices.append(all_vertices[i])
        
st_list = []
for s in req_vertices:
    for t in req_vertices:
        if s < t:
            st_list.append((s,t))

edges = list(edge_weights.keys())
for e in edges:
    edge_weights[((e[1],e[0]))] = edge_weights[e]


ng = nx.Graph()
f = open(dir + '/edges.txt', 'r')
for l in f:
    l = l.strip().split(' ')
    ng.add_edge(int(l[0]),int(l[1]))
    ng.add_edge(int(l[1]),int(l[0]))

H = steiner_primal_dual(ng, st_list, edge_weights)

nx.draw(H, with_labels = True)
plt.show()







    
