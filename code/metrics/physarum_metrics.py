import networkx as nx
import scipy
import random
import operator
import numpy as np

#m1: length distribution (aka, cost)

def length_distr(graph,underlying_graph = None):
	# the nodes are taken from the underlying graph
	# (this is useful in the case where we rectify the graph but we want to measure the lengths with using the rectified
	# edges) [graph_reduction][2nd_filtering]

	if underlying_graph==None:
		underlying_graph = graph

	list_of_lengths = []
	for edge in underlying_graph.edges():
		L=0
		v=edge[0]
		u=edge[1]
		path = nx.shortest_path(graph, v, u)
		#print(path)
		assert len(path)<3 #loops
		for i in range(len(path)-1):
			 w1 = path[i]
			 w2 = path[i+1]
			 p1 = graph.nodes[w1]['pos']
			 p2 = graph.nodes[w2]['pos']
			 l = scipy.spatial.distance.euclidean(p1,p2)
			 L=L+l
		list_of_lengths.append(L)
	#list_of_lengths = [d['length'] for u,v,d in graph.edges(data=True)]
	return list_of_lengths

#m2: degree distribution

def degree_distr(G):
	return [G.degree(node) for node in G.nodes()]

#m3: minimum distance

def minimum_distance(G,weight='tdens'):

	weight_list = []

	# computing all the shortest paths (weighted)
	asp = dict(nx.all_pairs_dijkstra(G,weight=weight))
	
	nodes = list(G.nodes())
	for i in range(len(nodes)-1):
		node1 = nodes[i]
		for j in range(i+1,len(nodes)):
			node2 = nodes[j]
			len_path_ = asp[node1][0][node2] #0 to access the 'length' of the path's dictionary 
			'''
			if len_path_ != 1:
				for i in range(len(path_) - 1):
					w1 = path_[i]
					w2 = path_[i + 1]
					p1 = G.nodes[w1]['pos']
					p2 = G.nodes[w2]['pos']
					l = scipy.spatial.distance.euclidean(p1, p2)
					L += l
			'''
			weight_list.append(len_path_)


	#output: minimum distance min len (path (u,v)) for all pairs u,v

	return weight_list

#m4: breaking_probabilty (~fault tolerance)


def efficiency(length_list):
	return np.array([1./l for l in length_list]).sum()/float(len(length_list))

def breaking_prob(G):


	edges = list(G.edges())

	ncc_diff_list = []

	for v, u in edges:
		ncc_G = len(list(nx.connected_components(G)))

		G_copy = G.copy()
		G_copy.remove_edge(v, u)
		ncc_G_copy = len(list(nx.connected_components(G_copy)))

		ncc_diff = ncc_G_copy - ncc_G  # either 1 or 0

		assert ncc_diff in [0, 1]

		ncc_diff_list.append(ncc_diff)

	#output: probability of disconnecting the graph by removing one edge

	return 1 - sum(ncc_diff_list) / len(ncc_diff_list)

#m5: transport efficiency

def transport_efficiency(G):
	_, length_list = minimum_distance(G)
	N = len(G.nodes())
	rec_lenghts = [1/elem for elem in length_list]
	te=sum(rec_lenghts)*(1/(N*(N-1)))
	return te

#m6: fault tolerance (targeted or random, node or edge)

def fault_tolerance(G, elem_type, removing_flag, removal_percentage, weight_flag ='weight',rseed=0, deg = None,  bn = None):
    
    rng = np.random.RandomState(seed=rseed)
    
    dg_flag = False
    
    if deg is not None:
        dg_flag = True
        #print('Not None')
    
    bn_flag = False
    if bn is not None:
        bn_flag = True

    G_copy = G.copy()
    #print('number of nodes',len(G_copy.nodes))

    # defining the 'bucket' (i.e. set of options to remove)
    if elem_type == 'node':
        bucket = list(G.nodes())
    elif elem_type == 'edge':
        bucket = list(G.edges())
    else:
        print('elem_type not defined!')

    bucket = np.array(bucket)
    N = len(bucket)

    # defining the number of removals (a % of the things in the bucket)

    if removal_percentage>1:

        rem_bound = removal_percentage

    else:

        rem_bound = int(removal_percentage * N)
        
    if rem_bound == 0: # for very small percentages
        
        rem_bound+=1

    No = len(G.nodes())

    if No-rem_bound !=0:

        # removing by cases:

        ## random: we take 'rem_bound' elements from the bucket

        if removing_flag == 'random':

            idx_to_be_removed = rng.choice(len(bucket), rem_bound)
            
            to_be_removed = bucket[idx_to_be_removed]

        ## targeted: we remove the 'most important' elements

        elif removing_flag == 'targeted':


            if elem_type == 'node': # we remove the nodes with highest degree
                
                if 'betweenness' in weight_flag or 'degree' in weight_flag:
                    
                    if not bn_flag:
                        bn = nx.betweenness_centrality(G) 
                    if not dg_flag:
                        deg = nx.degree_centrality(G)
                    
                    deg_sorted_rem_bound= {key: deg[key] for key in sorted(deg, key=deg.get, reverse=True)[:rem_bound]}
                    
                    bn_sorted_rem_bound = {key: bn[key] for key in sorted(bn, key=bn.get, reverse=True)[:rem_bound]}
                    
                    #print(deg_sorted_rem_bound)
                    
                    if weight_flag == 'degree':      
                        
                        to_be_removed = list(deg_sorted_rem_bound.keys())
                        
                        #print('different-value elements:',len(set(deg.values())),'/',len(deg.keys()))
                        
                    elif weight_flag =='degree+betweenness':

                        # getting the min degree of the candidates (why? so we dont get these rem_bound amount of nodes sorted just by labels)

                        min_deg = min(deg_sorted_rem_bound.values())

                        #print('len before adding the equal-degree ones',len(deg_sorted_rem_bound),'of',rem_bound)

                        equal_degree = {node: deg[node] for node in deg.keys() if deg[node]==min_deg}

                        top_sorted_without_equal_degree = [node for node in deg_sorted_rem_bound.keys() if node not in equal_degree.keys()]

                        equal_degree_sorted_by_bn =  [key for key in sorted(bn, key=bn.get, reverse=True)[:rem_bound - len(top_sorted_without_equal_degree)]] #the reamining rem_bound - len (taken) nodes

                        to_be_removed = top_sorted_without_equal_degree + equal_degree_sorted_by_bn

                    elif weight_flag == 'betweenness':
                        
                        to_be_removed = list(bn_sorted_rem_bound.keys())
                        
                        #print('different-value elements:',len(bn.values()),'/',len(bn.keys()))
                        
                    else:
                        
                        print('weight_flag not defined')
                        
                elif 'clustering' in weight_flag:

                    cl = nx.clustering(G)

                    cl_sorted_rem_bound= {key: cl[key] for key in sorted(cl, key=cl.get, reverse=True)[:rem_bound]}

                    to_be_removed = list(cl_sorted_rem_bound.keys()) 
                    
                    #print('different-value elements:',len(cl.values()),'/',len(cl.keys()))

                else:
                    print('weight_flag not defined')

                #to_be_removed = rem

            elif elem_type == 'edge': # we remove the edges with largest weight (eihter 'weight' or 'length')

                weights = {edge: G.edges[edge][weight_flag] for edge in bucket}

                weights_sorted_rem_bound = {key: weights[key] for key in sorted(weights, key=weights.get, reverse=True)[:rem_bound]}

                to_be_removed = list(weights_sorted_rem_bound.keys())


        else:

            print('removing_flag not defined!')

        #print('tbr',to_be_removed)


        if elem_type == 'node':

            G_copy.remove_nodes_from(to_be_removed)

        elif elem_type == 'edge':

            G_copy.remove_edges_from(to_be_removed)

        else:
            print('wrong input')
            
        G_copy.remove_nodes_from(list(nx.isolates(G_copy)))

        N1 = len(G_copy.nodes())
        
        try:
            largest_cc = max(nx.connected_components(G_copy), key=len) # this is a list of nodes
        except:
            print('graph is too small')
            largest_cc = []
            N1 = 1
        # if elem_type == 'edge': #using the same removal_percentage given for the edges causes inconsistencies

        # 	removal_percentage = (No-len(G_copy.nodes()))/No
        # 	print('No',No, 'removed',No-len(G_copy.nodes()), '%',removal_percentage)
        # 	print(removal_percentage)

        # denominator: the remaining number of NODES after removing

        # denom =int((1-removal_percentage)*No)
        denom = N1
        
        #print('number of elements removed:',len(to_be_removed))

        # debugging:

        cc = list(nx.connected_components(G_copy))
        len_cc = [len(c) for c in cc]

        #assert No >= sum(len_cc) + round(removal_percentage*No) - 1

        gcc = len(largest_cc) / denom
        
        #print('lcc',len(largest_cc),'den',denom)

    else:
        print('no nodes left.')
        gcc = 0

    #print('gcc',gcc)

    return gcc, G_copy


def rescaling_graph(G):
    
    edges = list(G.edges())
    
    try:
        
        weights = [G.edges[e]['tdens'] for e in edges]
        
        #print('works!')
        
        flag = 'tdens'
        
    except:

        weights = [G.edges[e]['weight'] for e in edges]

        flag = 'weight'
        
    maxw = max(weights)
    
    for e in edges:
        
        G.edges[e][flag] = G.edges[e][flag] / maxw
        
    assert max([G.edges[e][flag] for e in edges]) == 1
    
    return G
    
    
