import pickle as pkl
import decimal
import scipy
import itertools
from decimal import Decimal
import networkx as nx
import numpy as np

def loading_graph(folder_name, funct, t_string, graph_type,flag, min_string=None, btns_factor_string=None):
    if flag=='original':
        path_=folder_name+'/'+funct+'/'+funct+'_graph_t'+t_string+'_graph'+str(graph_type)+'.dat'
    elif flag=='simplification':
        path_=folder_name+'/'+funct+'/BP_simpl_'+funct+'_t'+t_string+'_graph'+str(graph_type)+'_min'+min_string+'_btns_factor'+btns_factor_string+'.dat'
    with open(path_, 'rb') as file:
        Graph=pkl.load(file)
    return Graph        

'''
def weight_part(G, M):
    
    This function computes the difference between the sum of the weights of a graph and the sum of the weights
    of the triangulation.
    :param G: Graph.
    :param M: Total weight of the partition.
    :return:
        diff_weight: difference between the weights.
    
    G_weight=sum([G.edges[edge[0],edge[1]]['weight'] for edge in G.edges(data=True)])
    diff_weight= abs(G_weight-M)/M
    return diff_weight
'''

def lengths_of_graph(G):
    '''
    This function computes the total length of a graph.
    :param G: Graph.
    :return:
        total_length: total length of the graph.
    '''
    length_dict={}
    i=1
    for edge in G.edges():
        length_dict[i]=scipy.spatial.distance.euclidean(G.nodes[edge[0]]['pos'],G.nodes[edge[1]]['pos'])
        i+=1
    total_length = sum([value for value in length_dict.values()]) 
    return total_length

def connect_part(G):
    '''
    This computes the number of "disconnections" of a graph.
    :param G: Graph.
    :return:
        output: number of connected components - 1.
    '''
    return len(list(nx.connected_components(G)))-1

def q_measure_old(G,totalW, number_of_edges_of_ugs, number_of_nodes_of_ugs):
    q_m=weight_part(G, totalW)+lengths_of_graph(G)/number_of_edges_of_ugs + connect_part(G)/number_of_nodes_of_ugs
    return q_m

def old_label(node_new,G_old,G_new):
    '''
    Given a node (with a label) in a graph, this function outputs the label of this node in another graph.
    :param node_new: node in G_new.
    :param G_old: labeled graph.
    :param G_new: labeled graph.
    :return:
        n: label of node in G_old.
        p: geometric location of the node (the same in both graphs).
    '''
    p =G_new.nodes[node_new]['pos']
    for n in G_old.nodes():
        if G_old.nodes[n]['pos'][0]==p[0] and G_old.nodes[n]['pos'][1]==p[1]:
            w = G_old.nodes[n]['weight']
            return n, p, w
        
def relabeling(G_new,G_old):
    '''
    Given two graphs G_new and G_old, this functions assign the labels of G_old, to the nodes of G_new using
    geometric locations as references.
    :param G_new: labeled graph.
    :param G_old: labeled graph.
    :return:
        G_new_relabeled: graph with new labels.
    '''
    test_edge = list(G_new.edges())[0] #{'w':--,'flux':--,...}
    try:
        w = G_new.edges[test_edge]['tdens']
        weight_label = 'tdens'
    except:
        weight_label = 'weight'

    G_new_relabeled=nx.Graph()
    for edge in G_new.edges():
        label1,p1, w1 = old_label(edge[0],G_old,G_new)
        label2,p2, w2 = old_label(edge[1],G_old,G_new)
        if weight_label =='weight':
            G_new_relabeled.add_edge(label1,label2, 
            weight=G_new.edges[edge[0],edge[1]][weight_label])
        else:    
            G_new_relabeled.add_edge(label1,label2, 
                tdens=G_new.edges[edge[0],edge[1]]['tdens'],
                #flux=G_new.edges[edge[0],edge[1]]['flux']
                )
        G_new_relabeled.nodes[label1]['pos']=p1
        G_new_relabeled.nodes[label2]['pos']=p2
        G_new_relabeled.nodes[label1]['weight']= w1
        G_new_relabeled.nodes[label2]['weight'] = w2
    return G_new_relabeled

def reweighting(G_new,G_old):
    '''
    This function assings the weights of the edges of a graph to the edges of another graph using
    geometric locations as references.
    :param G_new: labeled graph.
    :param G_old: labeled graph.
    :return:
        G_new: graph with new weights.
    '''
    for edge in G_new.edges():
        G_new.edges[edge[0],edge[1]]['tdens']=G_old.edges[edge[0],edge[1]]['tdens']
        G_new.edges[edge[0],edge[1]]['flux']=G_old.edges[edge[0],edge[1]]['flux']
    return G_new
    
def partition(N):
    '''
    Given an integer, this gives the set of coordinates of the (N-1)^2 - regular partition of the set [0,1]^2.
    :param N: integer.
    :return:
        x_ y_coord: list of coordinates.
        new_grid: graph whose nodes are the vertices of the elements of the partition and its edges their sides.
    '''
    grid=nx.grid_2d_graph(N,N)
    new_grid=nx.convert_node_labels_to_integers(grid)
    x_coord=list(np.linspace(0,1,N))
    y_coord=list(np.linspace(0,1,N))
    x_y_coord=list(itertools.product(x_coord,y_coord))
    for node in new_grid.nodes():
        new_grid.nodes[node]['pos']= x_y_coord[node]
    return x_y_coord, new_grid

'''
def partition_set(N):
    
    This returns a dict with some indexes for the elements of the partition. An element of 
    the partition is determined by its four 'corners'
    
    part_dict={}
    x_coord=list(np.linspace(0,1,N))
    y_coord=list(np.linspace(0,1,N))
    x_y_coord=list(itertools.product(x_coord,y_coord))
    N=np.sqrt(len(x_y_coord))
    number_of_elements_in_partition=(N-1)*(N-1)
    #print(number_of_elements_in_partition)
    for i in range(int(number_of_elements_in_partition)):
        #print(i)
        #print(x_y_coord[i])
        part_dict[i+1]=[x_y_coord[i],x_y_coord[i+1],x_y_coord[i+int(N)],x_y_coord[i+1+int(N)]]
    return part_dict
'''


def partition_set(N):
	'''
	Generates a regular partition of the set [0,1]^2.
	:param N: number of divisions of the interval [0,1] (endpoints count as divisions).
	:return:
		part_dict: dictionary, s.t., part_dict[key]=[(x1,y1),...,(x4,y4)].
		 dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2,n3 and n4 needs to be defined as dict_seq[key]=[n1,n2,n3,n4].
		node2box_index:  given a node, node2box_index[node] returns the indices of all the elements of the partition s.t.,
		node is a vertex of these elements.
	'''
	x_y_coord, _ = partition(N)
	part_dict = {}
	dict_seq = {}
	node2box_index = {}

	number_of_elements_in_partition = (N - 1) * (N - 1)
	for node in range(N ** 2):
		# print(node)
		node2box_index[node] = []
	# print('n',number_of_elements_in_partition)
	nbox = 1
	i = 0
	while nbox <= number_of_elements_in_partition:
		# print('congr',i % N)
		if i % N != (N - 1):
			# print(i)
			part_dict[nbox] = [x_y_coord[i], x_y_coord[i + 1], x_y_coord[i + int(N)], x_y_coord[i + 1 + int(N)]]
			dict_seq[nbox - 1] = [i, i + 1, i + int(N), i + 1 + int(N)]
			# print(dict_seq)
			for k in dict_seq[nbox - 1]:
				node2box_index[k].append(nbox - 1)
			nbox += 1
		i += 1
	return part_dict, dict_seq, node2box_index

def local_graph_weight_sum(G, list_of_coord, binary = False):
    """
    This script computes the local weight of a graph inside a member of the partition whose location is given by
    list_of_coord.
    Edges with both endpoints included in the element of the part. contribute totally to the overall weight
    Edges with exactly one endpoint included in the elem of the part. contribute half weight to the overall weight
    Edges with no endpoints in the elem. don't contribute to the sum.
    :param G: Graph.
    :param list_of_coord: coordinates of an element of the partition.
    :return:
        localW: local weight of G.
    """
    test_edge = list(G.edges())[0] #{'w':--,'flux':--,...}
    try:
        w = G.edges[test_edge]['tdens']
        weight_label = 'tdens'
    except:
        weight_label = 'weight'


    localW=0
    x1_cell=list_of_coord[0][0]
    x2_cell=list_of_coord[3][0]
    y1_cell=list_of_coord[0][1]
    y2_cell=list_of_coord[1][1]
    for edge in G.edges():
        v=edge[0]
        w=edge[1]
        vx=G.nodes[v]['pos'][0]
        vy=G.nodes[v]['pos'][1]
        wx=G.nodes[w]['pos'][0]
        wy=G.nodes[w]['pos'][1]
        bool_membership=[ x1_cell <= vx <= x2_cell and  y1_cell <= vy <= y2_cell,
                         x1_cell <= wx <= x2_cell and  y1_cell <= wy <= y2_cell ]
        if binary: 
            w_edge = 1
        else:
            w_edge = G.edges[v,w][weight_label]

        delta_edge=1/2*sum(bool_membership)
        localW=localW+delta_edge*w_edge
    return localW

def local_triang_weight_sum(G, list_of_coord, binary = False, threshold = 0):
    '''
    This script computes the amount of weight inside a member of the partition. 
    The weights are given by the barycenter graph.
    :param G: Graph.
    :param list_of_coord: coordinates of an element of the partition.
    :return:
        localTW: local weight of partition.
    '''
    localTW=0
    x1_cell=list_of_coord[0][0]
    x2_cell=list_of_coord[3][0]
    y1_cell=list_of_coord[0][1]
    y2_cell=list_of_coord[1][1]
    for node in G.nodes():
        vx=G.nodes[node]['pos'][0]
        vy=G.nodes[node]['pos'][1]
        if x1_cell<= vx <= x2_cell and y1_cell <= vy <= y2_cell:
            w_node=G.nodes[node]['weight']
            if binary: 
                if w_node>threshold:
                    w_node = 1
                else:
                    w_node = 0
            else:
                pass
            #print(w_node)
            localTW=localTW+ w_node
    return localTW

def q_measure(G,opt_triang_G, part_dict, flag_dist,min_, binary = False):
    '''
    This script computes the quality measure of a graph G built from a optimal triangulation w.r.t., a partition.
    :param G: Graph.
    :param opt_triang_G: Graph representing the optimal triangulation.
    :param part_dict: dictionary, s.t., part_dict[key]=[(x1,y1),...,(x4,y4)].
    :param flag_dist: 'l2' for l2-distance; 'l1' for l1-distance.
    :return:
        q_measure: value of the measure.
        q_weight: value of the weight-related part.
        L: length of the Graph.
        triang_weight_vector: vector containing the value of the weights of the optimal triangulation inside each
        element of the partition.
        weight_difference_vector: vector containing the value of the quality measure inside each element of the
        partition.
    '''
    max_=max([opt_triang_G.nodes[node]['weight'] for node in opt_triang_G.nodes()])
    '''
    # Removing non-useful weights from triangulation
    for node in opt_triang_G.nodes:
        w = opt_triang_G.nodes[node]['weight']
        if w <= min_ * max_:
            opt_triang_G.nodes[node]['weight']=0
    '''
    weight_difference_vector=[]
    triang_weight_vector=[]
    for coordinates in part_dict.values():
        #computing the weight of the opt triang in the coordinates elem. of the part
        w_triang = local_triang_weight_sum(opt_triang_G, coordinates, binary, min_*max_)
        #computing the local weight of the G in the coordinates elem. of the part
        w_graph = local_graph_weight_sum(G, coordinates, binary)
        #appending to the lists
        triang_weight_vector.append(w_triang)
        weight_difference_vector.append(abs(w_triang-w_graph))
    L = lengths_of_graph(G)
    part_factor = 1/len(part_dict.keys())
    if flag_dist == 'l2':    
        q_weight = part_factor*np.linalg.norm(weight_difference_vector)
        q_measure = q_weight + L
    elif flag_dist == 'l1':
        q_weight = part_factor*sum([abs(x) for x in weight_difference_vector])
        q_measure =q_weight+L

    return q_measure, q_weight, L, triang_weight_vector, weight_difference_vector



ignore=True

if not ignore:

    def loading_bar_structure(folder_name, funct,t):
        '''
        Importing the graph structure
        We take only the bareycenter and then we compute the edges.
        '''
        graph_coord_triang,graph_coordinates,n_nodes=get_baryc(folder_name)
        '''
        Saving the positions of the barycenters in a dict(to plot!)
        '''
        bar_pos=bar2dict(graph_coordinates,n_nodes)
        '''
        Extracting the weights (opt_tdens, opt_pot, opt_flux)
        '''
        if funct=='tdens':
            opt_tdens=extracting_weights(folder_name,'output/result/opt_tdens.dat')
        if funct=='flux':
            opt_tdens=extracting_weights(folder_name,'output/result/opt_tdens.dat')
            opt_pot=extracting_weights(folder_name,"output/result/opt_nrm_grad_avg.dat")
        '''
        Saving the tdens (or pot or flux) in a dict and storing the max
        '''
        if funct=='tdens':
            dict_weights_tdens,weights_tdens=weight2dict(opt_tdens,'output/result/opt_tdens.dat')
        if funct=='flux':
            dict_weights_tdens,weights_tdens=weight2dict(opt_tdens,'output/result/opt_tdens.dat')
            dict_weights_pot,weights_pot=weight2dict(opt_pot,"output/result/opt_nrm_grad_avg.dat")
            dict_weights,weights=prod_dict(dict_weights_tdens,dict_weights_pot,weights_tdens,weights_pot)
        else:
            dict_weights=dict_weights_tdens
            weights=weights_tdens
        max_=max(weights)
        '''
        Completing with zeros
        '''
        dict_weights=completing_with_zeros(dict_weights,bar_pos)
        '''
        Defining the graph with the attributes (pos and tdens)
        '''
        X_func=nx.Graph()
        #X_func.add_nodes_from(bar_pos.keys())
        #print(X_func.nodes())
        for n, p in bar_pos.items():
            w=dict_weights[int(n)]
            if w>= t*max_:
                X_func.add_node(n)
                X_func.node[n]['pos'] = p
                X_func.nodes[n]['weight'] = dict_weights[int(n)]


        pos=nx.get_node_attributes(X_func,'pos')

        return X_func







