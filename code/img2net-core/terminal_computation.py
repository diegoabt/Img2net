import time
import networkx as nx
import numpy as np
import operator
import itertools
#--------------------
import quality_measure
#--------------------

def pixel_circle_counter(center,R, partition_dict, G):
    intersecting_circle=[]
    for key in G.nodes():
        corners = partition_dict[key+1]
        bool_list = [(center[0]-R<=point[0]<=center[0]+R) & (center[1]-R<=point[1]<=center[1]+R) for point in corners]
        if 0<sum(bool_list)<4:
            intersecting_circle.append(key+1)
    return intersecting_circle


def single_box_finder(intersection_list, partition_dict):  # under construction
    # placing the first square in the first box
    box = [intersection_list[0]]
    box_coord = partition_dict[intersection_list[0]]
    # print('bc',box_coord)
    for key in intersection_list[1:]:
        # testing if there's an intersection
        # print(key)
        chosen_coord = partition_dict[key]
        # print('cc',chosen_coord)
        A = len(box_coord)
        # print('a',A)
        # print('--u',chosen_coord+box_coord)
        # print('--su',set(chosen_coord+box_coord))
        AunionB = len(set(chosen_coord + box_coord))
        # print('ab',AunionB)
        if AunionB != A + 4:
            box_coord = list(set(box_coord + chosen_coord))
            box = box + [key]
    return box, box_coord

def cluster_counter(intersection_list,partition_dict):
    box={}
    len_intersect = len(intersection_list)
    box_total=0
    while len_intersect>0:
        box,_=single_box_finder(intersection_list,partition_dict)
        intersection_list = [n for n in intersection_list if n not in box]
        len_intersect=len(intersection_list)
        box_total+=1
    return box_total

def branch_counter(intersection_list,partition_dict):

    box_total = cluster_counter(intersection_list,partition_dict)
    
    if box_total==2 or box_total==0:
        #path
        nbr = 0
    elif box_total==1:
        #endpoint
        nbr=-1
    else:
        #bifurcation,trifurcations,...
        nbr=1
    return nbr


def terminal_finder(D, partition_dict, G):
    
    print('graph size:' ,len(G.nodes) ,' nodes' ,len(G.edges) ,' edges.')
    start_time = time.time()
    terminal_list =[]
    color_nbr =[]
    nodes_for_correction = []

    # bif graph, i.e., terminal mapping

    nbr_graph =nx.Graph()

    # betweenness centrality

    bn =nx.betweenness_centrality(G)

    # time

    print("bn centralities: --- %s second(s) ---" % (time.time() - start_time))
    
    start_time = time.time()

    # closeness centrality

    cn =nx.closeness_centrality(G)

    # time

    print("cn centralities: --- %s second(s) ---" % (time.time() - start_time))

    start_time = time.time()

    # building the centers of the masks

    number_of_masks =round( 1 /(D))

    print('number of masks:' , 1 /D ,number_of_masks)

    
    x_coord_mask =list(np.linspace( D /2 , 1 - D /2 ,number_of_masks))
    x_y_coord_mask =list(itertools.product(x_coord_mask ,x_coord_mask))

    # sweeping the masks

    N=0
    for center in x_y_coord_mask:

        correction_cluster= []
        N+=1

        # get the indices of the intersected pixels

        intersection_list = pixel_circle_counter(center,D/2, partition_dict, G)

        # get the *bif*

        nbr = branch_counter(intersection_list,partition_dict)

        # get all the nodes in the mask

        nodes_in_the_mask=[node for node in G.nodes()
                           if (
                             center[0]-D/2<= G .nodes[node]['pos'][0]<center[0]+D/2)
                           & (center[1]-D/2<= G .nodes[node]['pos'][1]<center[1]+D/2)
                           ]
        
        # add them to the terminal mapping graph

        nbr_graph.add_nodes_from(nodes_in_the_mask)

        # compute the subgraph G_C

        subGraph=G.subgraph(nodes_in_the_mask)

        # Case on *bif*

        if nbr==-1 and len(nodes_in_the_mask)>0:

            # the least closeness centrality in the mask

            cn_in_the_mask = {n:cn[n] for n in nodes_in_the_mask}

            term = min(cn_in_the_mask.items(), key=operator.itemgetter(1))[0]

            # add to the list of terminals

            terminal_list.append(term)

        elif nbr>0 and len(nodes_in_the_mask)>0:

            # get the maximum betweenness centrality in the filter

            cn_in_the_mask = {n:bn[n] for n in nodes_in_the_mask}

            term = max(cn_in_the_mask.items(), key=operator.itemgetter(1))[0]
            
            # NOTE: this is not appended, maybe because we are using the conjecture

            # terminal_list.append(term)

            # compute the closeness center for each one

        # add this node to the terminal l is t

        elif nbr==0:

        	# get connected components of G_c
            
            ccom = list(nx.connected_components(subGraph))

            lcc = len(ccom )

            if lcc==1:
                # print('nothing special in this 0-branch')
                pass
            elif lcc>1:

            	# very likely to be always 2, if not 1.

                #nodes_for_correction[N]=[ ]

                # increase the size of the mask 10%

                Dn=D*(1.1)

                # get the new set of nodes

                nodes_in_the_extended_mask=[node for node in G.nodes()
                                            if (center[0]-Dn/2<=G. nodes[node]['pos'][0]<center[0]+Dn/2)
                                            & (center[1]-Dn/2<=G. nodes[node]['pos'][1]<center[1]+Dn/2)
                                            ]

                # get the new G_C for this aumented mask

                subGraph_extended =G.subgraph(nodes_in_the_extended_mask)

                # get the connected components

                ccom_extended = list(nx.connected_components(subGraph_extended))
                
                lcc_e = len(ccom_extended)

                # if they have the same length, then ... 
                if True:  # lcc==lcc_e:

                	# add to the set of terminals the nodes with the minimum betweenness centrality

                    for elem in ccom:

                        elem_subgraph = subGraph.subgraph(elem)

                        bn_sub = {n:bn[n] for n in elem_subgraph.nodes()}

                        term = min(bn_sub.items(), key=operator.itemgetter(1))[0]

                        terminal_list.append(term)

                        # add the nodes to the set \mathcal{E} for edge correction

                        correction_cluster.append(term)
                        
                    nodes_for_correction.append(correction_cluster)

        # coloring

        for node in nodes_in_the_mask:

            nbr_graph.nodes[node]['pos'] = G.nodes[node]['pos']

            #add the color depending on its nbr

            if nbr == -1: #end-point
                color_nbr.append('red')
            elif nbr == 0:#path
                color_nbr.append('blue')
            else:#bifurcation
                color_nbr.append('green')

    print("rest: --- %s second(s) ---" % (time.time() - start_time))
    return nbr_graph, color_nbr, terminal_list, nodes_for_correction, number_of_masks