import os
import pickle as pkl
import networkx as nx
import scipy
from shutil import copyfile
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import numpy as np
#----------------------------------------------------------
import quality_measure
import source_sink_generator
import terminal_computation
import pre_extraction
import utils
#---------------------------------------------------------


def concatenate(lists):
    '''
    This concatenates all the lists contained in a list.
    :param lists: a list of lists.
    :return:
        new_list: concatenated list.
    '''
    new_list = []
    for i in lists:
        new_list.extend(i)
    return new_list


def terminals_from_cont(Graph, source_flag, sink_flag, btns_factor_source, btns_factor_sink, terminal_criterion='branch_convex_hull+btns_centr'):
    '''
    Computation of source and sink nodes. This script uses information about the inputs of the DMK solver.
    There are three criteria for the selection of the nodes.
    :param Graph: a networkx graph to be filtered.
    :param source_flag: flag used to define the source region in the continuous problem.
    :param sink_flag: flag used to define the sink region in the continuous problem.
    :param btns_factor_source: threshold for the nodes in the source region (see more in paper).
    :param btns_factor_sink: threshold for the nodes in the sink region (see more in paper).
    :param terminal_criterion: 'branch_convex_hull+btns_centr' (combination of btns centr and convex hull by branches),
    'whole_convex_hull+btns_centr' (combination of btns centr and convex hull of the source and sink regions),
    'btns_centr' (only btns centr).

    :return:
        possible_terminals_source: for each i, possible_terminals_source[i]= "sources" of i-th cc.
        possible_terminals_sink: for each i, possible_terminals_sink[i]= "sources" of i-th cc.
    '''


    bn = nx.betweenness_centrality(Graph, normalized=True)

    # Defining source-sink nodes of the graph (nodes inside the source or sink of the continuous DMK)
    if source_flag == 'rect_cnst':
        nodes_in_source = [node for node in Graph.nodes() if
                           source_sink_generator.source_rect_cnst_test(Graph.nodes[node]['pos'][0],
                                  Graph.nodes[node]['pos'][1])]
    else:
        nodes_in_source = [node for node in Graph.nodes() if
                           source_sink_generator.fsource(Graph.nodes[node]['pos'][0], Graph.nodes[node]['pos'][1],
                                   source_flag) != 0]
    if sink_flag == 'rect_cnst':
        nodes_in_sink = [node for node in Graph.nodes() if
                         source_sink_generator.sink_rect_cnst_test(Graph.nodes[node]['pos'][0],
                              Graph.nodes[node]['pos'][1])]
    else:
        nodes_in_sink = [node for node in Graph.nodes() if
                         source_sink_generator.fsink(Graph.nodes[node]['pos'][0], Graph.nodes[node]['pos'][1],
                               sink_flag) != 0]

    # min bn inside the source and sink

    max_bn_source = max([bn[node] for node in nodes_in_source])
    max_bn_sink = max([bn[node] for node in nodes_in_sink])

    # Defining source-sink candidates based only on btns

    kind_of_leaf_nodes_source = [key for key in nodes_in_source if
                                 bn[key] <= btns_factor_source * max_bn_source]  # *(min_bn_source)+.0001]
    kind_of_leaf_nodes_sink = [key for key in nodes_in_sink if
                               bn[key] <= btns_factor_sink * max_bn_sink]  # *(min_bn_sink+.0001)]

    # Defining the subgraphs induced by the candidates

    sub_source = Graph.subgraph(kind_of_leaf_nodes_source)
    sub_sink = Graph.subgraph(kind_of_leaf_nodes_sink)

    # Removing repeated labels

    possible_terminals_source = set(kind_of_leaf_nodes_source)
    possible_terminals_sink = set(kind_of_leaf_nodes_sink)


    if terminal_criterion != 'btns_centr':
        # Getting the coordinates to compute convex hulls

        coordinates_in_source = np.asarray(
            [[Graph.nodes[node]['pos'][0], Graph.nodes[node]['pos'][1]] for node in
             nodes_in_source])
        coordinates_in_source_list = concatenate(list(coordinates_in_source))

        coordinates_in_sink = np.asarray(
            [[Graph.nodes[node]['pos'][0], Graph.nodes[node]['pos'][1]] for node in
             nodes_in_sink])
        coordinates_in_sink_list = concatenate(list(coordinates_in_sink))

        # If the number of coordinates (or nodes) is not more than 3, then no convex hull computation

        if len(coordinates_in_source) >= 4 and len(coordinates_in_sink) >= 4:

            # Computing convex hull for the branches

            if terminal_criterion == 'branch_convex_hull+btns_centr':
                source_hull = ConvexHull(coordinates_in_source)
                index_source_hull = np.asarray(source_hull.vertices)
                nodes_source_hull = []
                coord_source_hull = np.asarray(
                    [coordinates_in_source[node] for node in index_source_hull])
                for node in nodes_in_source:
                    x, y = Graph.nodes[node]['pos']
                    if x in coord_source_hull[:, 0] and y in coord_source_hull[:, 1]:
                        nodes_source_hull.append(node)

                sink_hull = ConvexHull(coordinates_in_sink)
                index_sink_hull = np.asarray(sink_hull.vertices)
                nodes_sink_hull = []
                coord_sink_hull = np.asarray([coordinates_in_sink[node] for node in index_sink_hull])
                for node in nodes_in_sink:
                    x, y = Graph.nodes[node]['pos']
                    if x in coord_sink_hull[:, 0] and y in coord_sink_hull[:, 1]:
                        nodes_sink_hull.append(node)
                possible_terminals_source = set(kind_of_leaf_nodes_source + nodes_source_hull)
                possible_terminals_sink = set(kind_of_leaf_nodes_sink + nodes_sink_hull)

            # Computing convex hull for all the nodes defined as candidates

            elif terminal_criterion == 'whole_convex_hull+btns_centr':  # not working!
                single_source_hull = ConvexHull(coordinates_in_source_list)
                single_index_source_hull = np.asarray(single_source_hull.vertices)
                nodes_source_hull = []
                single_coord_source_hull = np.asarray(
                    [coordinates_in_source_list[node] for node in single_index_source_hull])
                for node in nodes_in_source:
                    x, y = Graph.nodes[node]['pos']
                    if x in single_coord_source_hull[:, 0] and y in single_coord_source_hull[:, 1]:
                        nodes_source_hull.append(node)
                single_sink_hull = ConvexHull(coordinates_in_sink_list)
                single_index_sink_hull = np.asarray(single_sink_hull.vertices)
                nodes_sink_hull = []
                single_coord_sink_hull = np.asarray(
                    [coordinates_in_sink_list[node] for node in single_index_sink_hull])
                for node in nodes_in_sink:
                    x, y = Graph.nodes[node]['pos']
                    if x in single_coord_sink_hull[:, 0] and y in single_coord_sink_hull[:, 1]:
                        nodes_sink_hull.append(node)
                possible_terminals_source = set(kind_of_leaf_nodes_source + nodes_source_hull)
                possible_terminals_sink = set(kind_of_leaf_nodes_sink + nodes_sink_hull)


    return possible_terminals_source, possible_terminals_sink




def bifurcation_paths(G, terminals=[], weights_list = ['weight']):
    '''
    This script takes a filtered graph and reduces its paths (sequences of nodes with degree 2) to a single edge.

    :param G:  filtered graph (networkx graph).
    :param terminals: union of source and sink nodes.
    :return:
        G: reduced graph.
    '''

    G = G.copy()
    N = len(G.nodes())
    deg_norm = nx.degree_centrality(G)
    deg = {}
    for key in deg_norm.keys():
        deg[key] = int(round((N - 1) * deg_norm[key]))
    deg_3 = [node for node in G.nodes() if deg[node] >= 3 or deg[node] == 1 or node in terminals]

    G_wo_bifuc = G.copy()
    for node in deg_3:
        G_wo_bifuc.remove_node(node)
    cc = list(nx.connected_components(G_wo_bifuc))


    for comp in cc:
        comp = set(comp)
        neighs = {neigh for node in comp for neigh in G.neighbors(node) if neigh not in comp}
        assert len(neighs) <= 2
        # adding the weights
        nodes = list(neighs)+list(comp)
        
        
        subgraph = G.subgraph(nodes)
        
        subgraph = subgraph.copy()
        
        G.remove_nodes_from(comp)
        
        #print('sg',subgraph.edges(data=True))
        '''
        if len(neighs) == 1:
            v = list(neighs)[0]
            neighs = (v,v)

        G.remove_nodes_from(comp)

        if len(neighs) == 2: 
            neighs = tuple(neighs)
            G.add_edge(neighs[0],neighs[1])
        '''
        for weight in weights_list:
            #print(weight)
            #print('sg',subgraph.edges())
            w = 0
            l=0

            for edge in subgraph.edges:
                #print(edge,subgraph.edges[edge])

                w+=subgraph.edges[edge][weight]
                w1 = edge[0]
                w2 = edge[1]
                p1 = subgraph.nodes[w1]['pos']
                p2 = subgraph.nodes[w2]['pos']
                l+=scipy.spatial.distance.euclidean(p1, p2)
            
            neighs = tuple(neighs)
            
            if len(neighs) == 1:
                v = neighs[0]
                neighs = (v,v)

            

            if neighs not in G.edges():
                G.add_edge(neighs[0],neighs[1])
            
            G.edges[neighs]['length'] = l

            if 'weight' == weight:
                G.edges[neighs]['weight'] = w
            elif 'flux' == weight:
                G.edges[neighs]['flux'] = w
            elif 'tdens' == weight:
                G.edges[neighs]['tdens'] = w
            elif 'time' == weight:
                G.edges[neighs]['time'] = w
            #print(weight,'===',w)

    return G


def BP_solver(folder_name, index):
    '''
    This script executes the BP_solver (a muffe_sparse_opt/dmk_folder.py sequence)
    :param folder_name: folder path where the outputs will be stored. It should be written "./runs/folder_name".
    :param index: integer representing the connected component to which this algorithm is going to be applied.
    :return:
    '''

    # Generating the output folders
    try:
        os.mkdir('../otp_utilities/muffe_sparse_optimization/simplifications/' + folder_name + "/component" + str(
            index) + '/')
    except OSError:
        pass
        #print("Creation of the directory %s failed." % (folder_name + "/component" + str(index) + '/'))


    for folder in ["output", "output/result", "output/vtk", "output/linear_sys", "output/timefun"]:
        try:
            os.mkdir('../otp_utilities/muffe_sparse_optimization/simplifications/' + folder_name + "/component" + str(
                index) + '/' + folder)
        except OSError:
            pass
            #print("Creation of the directory %s failed." % (folder_name + "/component" + str(index) + '/' + folder))
    os.system(
        'cp ../otp_utilities/muffe_sparse_optimization/simplifications/muffa.fnames  ' + '../otp_utilities/muffe_sparse_optimization/simplifications/' + folder_name + "/component" + str(
            index) + '/')

    # Copying the par_files
    for file in ['decay', 'pflux', 'pmass']:
        os.system(
            'cp  ../otp_utilities/muffe_sparse_optimization/simplifications/par_files/' + file + '.dat ' + '../otp_utilities/muffe_sparse_optimization/simplifications/' + folder_name + "/component" + str(
                index) + '/input')

    # Executing the dmk_folder.py run
    continuous_path = os.path.abspath("./")
    #print('path in BP_solver', os.getcwd(), continuous_path)
    discrete_path = "../otp_utilities/muffe_sparse_optimization/simplifications/"


    os.chdir(discrete_path)
    #print('discrete', os.getcwd())
    command = "./dmk_folder.py run  " + folder_name[2:] + "/component" + str(
        index) + "  " + " muffa.ctrl  > outputs_dmk_d.txt"

    #print(command)
    os.system(command)
    os.chdir(continuous_path)





def filtering_from_image(small_G_filtered, beta_d, terminals, color_dict, partition_dict, weighting_method_simplification='ER',entries=[0],folder_name='./runs/no_name/',plotting=True):
    '''
    This takes as input a pre-extracted graph (obtained from an image) and filters it using filtering().
    :param small_G_filtered: pre-extracted graph.
    :param beta_d: beta input for the DMK solver.
    :param terminals: union of source and sink nodes.
    :param color_dict: dictionary s.t., color_dict[key]= real value for the key-th pixel.
        key is the index for the pixels in the resized image.
    :param partition_dict: dictionary, s.t., part_dict[key]=[(x1,y1),...,(x4,y4)].
    :param weighting_method_simplification: 'ER', 'IBP', 'BPW'.00
    :param entry: node index. terminals[entry] is the unique source node.
    :param folder_name: folder path where the outputs will be stored. It should be written "./runs/folder_name".
    :return:
        G_filtered: filtered graph (networkx graph).
    '''

    Graph = small_G_filtered.copy()
    BP_weights = 'BPtdens'
    min_ = .0001
    terminal_info =[terminals, entries]
    weight_flag='length'
    input_flag = 'image'

    #print(folder_name)
    folder_name='./runs/'+folder_name.split('/')[-1]
    print(folder_name)
    try:
        #print('Creating folder',folder_name.split('/')[-1])
        os.mkdir('../otp_utilities/muffe_sparse_optimization/simplifications/runs/'+folder_name.split('/')[-1])
    except:
        pass

    try:
        #print('Creating folder',folder_name)
        os.mkdir(folder_name)
    except:
        pass


    G_filtered, newGraph, ncc, possible_terminals_source, possible_terminals_sink, mapping, conv_report = filtering(Graph,
              beta_d,
                min_,
                folder_name+"/",
                terminal_info,
                weighting_method_simplification,
                BP_weights,
                weight_flag,
                input_flag)

    #print('filtering_from_image',G_filtered.edges(data=True))
    # Plotting
    ## Filtering on image
    if plotting:
        fig, ax = plt.subplots(1, 1, figsize=(15, 15))
        patches = []
        for key in partition_dict:
            square_edges = np.asarray([partition_dict[key][0]] + [partition_dict[key][2]] + [partition_dict[key][3]] + [
                partition_dict[key][1]] + [partition_dict[key][0]])
            s1 = Polygon(square_edges)
            patches.append(s1)
        p = PatchCollection(patches, alpha=.9, cmap='Wistia', linewidth=.0, edgecolor='b')

        colors = np.array(list(color_dict.values()))
        p.set_array(colors)
        ax.add_collection(p)

        plt.savefig(folder_name + '/image', transparent=False)

        fig, ax = plt.subplots(1, 1, figsize=(15, 15))
        patches = []
        for key in partition_dict:
            square_edges = np.asarray([partition_dict[key][0]] + [partition_dict[key][2]] + [partition_dict[key][3]] + [
                partition_dict[key][1]] + [partition_dict[key][0]])
            # print(square_edges)
            # print(len(square_edges))
            s1 = Polygon(square_edges)
            patches.append(s1)
        p = PatchCollection(patches, alpha=.9, cmap='Wistia', linewidth=.1, edgecolor='b')

        colors = np.array(list(color_dict.values()))
        p.set_array(colors)
        ax.add_collection(p)

        # plt.savefig('./component1/image', transparent=False)


        pos = nx.get_node_attributes(small_G_filtered, 'pos')
        nx.draw_networkx(small_G_filtered, pos, node_size=5, width=2, with_labels=False, edge_color='gray', alpha=0.5,
                         node_color='black', ax=ax)

        pos = nx.get_node_attributes(G_filtered, 'pos')
        nx.draw_networkx(G_filtered, pos, node_size=5, width=3, with_labels=False, edge_color='black', alpha=1,
                         node_color='black', ax=ax)

        for node in terminals:
            if terminals.index(node) in entries:
                color = 'green'
            else:
                color = 'red'
            x = small_G_filtered.nodes[node]['pos'][0]
            y = small_G_filtered.nodes[node]['pos'][1]
            circle1 = plt.Circle((x, y), .015, color=color, fill=False, lw=4)
            ax.add_artist(circle1)

        plt.savefig(folder_name + '/filtered_graph+image', transparent=False)

        ## Filtering with highlighted terminals

        fig, ax = plt.subplots(1, 1, figsize=(15, 15))

        pos = nx.get_node_attributes(G_filtered, 'pos')
        nx.draw_networkx(G_filtered, pos, node_size=5, width=3, with_labels=False, edge_color='black', alpha=1,
                         node_color='black', ax=ax)

        plt.savefig(folder_name + '/filtered_graph', transparent=False)

        fig, ax = plt.subplots(1, 1, figsize=(15, 15))

        pos = nx.get_node_attributes(G_filtered, 'pos')
        nx.draw_networkx(G_filtered, pos, node_size=5, width=3, with_labels=False, edge_color='black', alpha=1,
                         node_color='black', ax=ax)

        for node in terminals:
            if terminals.index(node) in entries:
                color = 'green'
            else:
                color = 'red'
            x = small_G_filtered.nodes[node]['pos'][0]
            y = small_G_filtered.nodes[node]['pos'][1]
            circle1 = plt.Circle((x, y), .015, color=color, fill=False, lw=4)
            ax.add_artist(circle1)

        plt.savefig(folder_name + '/filtered_graph+terminals', transparent=False)

    with open(folder_name + '/filtered_graph.pkl', 'wb') as file:
        pkl.dump(G_filtered, file)
    plt.close('all')

    return G_filtered, conv_report



def img_pre_extr2filtering(image_path, filter_size, beta_d, weighting_method_simplification='ER', entries=None, terminal_list = None):
    '''
    This takes as input a path containing the bfs approximation of the pre-extracted graph. This pre-extracted graph has
     been obtained from an image. The output is the filtered graph.
    :param image_path: string.
    :param filter_size: radious of the filters for the terminal selection process.
    :param weighting_method_simplification: 'ER', 'IBP', 'BPW'.
    :param beta_d: beta input for the DMK solver.
    :return:
        G_filtered: filtered graph (networkx graph).
    '''

    ## -- getting pre_extraction and bfs tree approximation

    if entries is None:
        entries = [0]

    if terminal_list is None:
        terminal_list = []

    last_word = image_path.split('/')[-1]
    new_folder_name = last_word.split('.')[0]
    saving_path = './runs/' + new_folder_name

    file = '/bfs_extracted_graph.pkl'

    with open(saving_path + file, 'rb') as file:
        G_bfs = pkl.load(file)

    file = '/real_image.pkl'

    with open(saving_path + file, 'rb') as file:
        color_dict = pkl.load(file)

    file = '/real_part_dict.pkl'

    with open(saving_path + file, 'rb') as file:
        partition_dict = pkl.load(file)

    file = '/extracted_graph.pkl'
    with open(saving_path + file, 'rb') as file:
        G_pre_extracted = pkl.load(file)

    ## -- removing leaves

    deletion_number = 1 #<----------------------------- !

    G_bfs_less_leaves = leaves_removal(G_bfs, deletion_number)

    ## -- getting edges for correction

    gbfs_threshold = 30 #<--------------------------- !
    gpe_threshold = 5 #<---------------------------- !

    edges, source_target_pair = edge_correction_v2(G_bfs_less_leaves, G_pre_extracted, gbfs_threshold,
                                                       gpe_threshold)
    nodes_for_correction = []

    for e in source_target_pair:
        nodes_for_correction.append(e[0])
        nodes_for_correction.append(e[1])

    ## -- getting terminals

    if terminal_list == None:

        nbr_graph, color_nbr, terminal_list, _, filter_number = terminal_computation.terminal_finder(
            filter_size,
            partition_dict,
            G_bfs)

        plotting = True

    else:

        print('terminals provided by user.')

        plotting = False

    ## -- getting the filtering

    terminal_list = list(set(terminal_list + nodes_for_correction))

    G_filtered = filtering_from_image(G_bfs, beta_d, terminal_list, color_dict, partition_dict,
                                      weighting_method_simplification,
                                      entries, saving_path)

    ## -- adding the missing edges

    G_filtered = quality_measure.relabeling(G_filtered, G_bfs)

    G_filtered.add_edges_from(edges)

    # fixing potential missing locations

    for node in G_filtered.nodes():
        G_filtered.nodes[node]['pos'] = G_bfs.nodes[node]['pos']

    ##### plotting terminal map ######

    if plotting:

        fig, ax = plt.subplots(1, 1, figsize=(15, 15))

        partition_dict_mask, _, _ = quality_measure.partition_set(filter_number + 1)
        patches = []
        for key in partition_dict_mask:
            square_edges = np.asarray(
                [partition_dict_mask[key][0]] + [partition_dict_mask[key][2]] + [partition_dict_mask[key][3]] + [
                    partition_dict_mask[key][1]] + [partition_dict_mask[key][0]])
            # print(square_edges)
            # print(len(square_edges))
            s1 = Polygon(square_edges)
            patches.append(s1)
        p = PatchCollection(patches, alpha=.5, linewidth=1, edgecolor='b', facecolor='white')
        ax.add_collection(p)

        patches2 = []
        for key in partition_dict:
            square_edges = np.asarray([partition_dict[key][0]] + [partition_dict[key][2]] + [partition_dict[key][3]] + [
                partition_dict[key][1]] + [partition_dict[key][0]])
            # print(square_edges)
            # print(len(square_edges))
            s1 = Polygon(square_edges)
            patches2.append(s1)
        p2 = PatchCollection(patches2, alpha=.4, cmap='YlOrRd', linewidth=.1, edgecolor='b')

        colors = np.array(list(color_dict.values()))
        p2.set_array(colors)
        ax.add_collection(p2)

        pos = nx.get_node_attributes(nbr_graph, 'pos')
        nx.draw_networkx(nbr_graph, pos, node_size=15, width=1.5, with_labels=False, edge_color='red', alpha=1,
                         node_color=color_nbr, ax=ax)

        pos = nx.get_node_attributes(G_bfs, 'pos')
        nx.draw_networkx(G_bfs, pos, node_size=0, width=1.5, with_labels=False, edge_color='black',
                         alpha=0.5, ax=ax)

        for i in range(len(terminal_list)):
            node = terminal_list[i]
            # color=color_nbr[i]
            if node in nodes_for_correction:
                color = 'green'
                size = 1
            else:
                color = 'black'
                size = .5
            x = G_bfs.nodes[node]['pos'][0]
            y = G_bfs.nodes[node]['pos'][1]
            circle1 = plt.Circle((x, y), .01, color=color, fill=False, lw=size)
            ax.add_artist(circle1)
            # label = ax.annotate(str(node), xy=(x, y), fontsize=12, ha="center")

        plt.savefig(saving_path + '/terminal_map.png')



    return G_filtered

def leaves_removal(G_bfs, deletion_number):

    G_bfs_less_leaves = G_bfs.copy()

    for deletion in range(deletion_number):

        # getting the leaves

        deg = nx.degree_centrality(G_bfs_less_leaves)

        N = len(G_bfs_less_leaves.nodes)

        for node in deg.keys():
            deg[node] = round(deg[node] * (N - 1))

        nodes_to_remove = [node for node in G_bfs_less_leaves.nodes() if deg[node] == 1]

        # removing them

        G_bfs_less_leaves.remove_nodes_from(nodes_to_remove)

    return G_bfs_less_leaves


def edge_correction_v1(G_bfs, G_pre_extracted,nodes_for_correction):

    edges = []

    for index in range(len(nodes_for_correction)):
        #
        nodes_list = nodes_for_correction[index]

        for i in range(len(nodes_list) - 1):

            source = nodes_list[i]
            target = nodes_list[i + 1]

            len_ = nx.shortest_path_length(G_bfs, source=source, target=target)

            if len_ > 20:
                # print('=',source,target,'l',len_)
                # print('1.something to add!')
                nodes_to_add = nx.shortest_path(G_pre_extracted,
                                                source=source,
                                                target=target)
                if len(nodes_to_add) < 15:
                    # N+=1
                    # print('2.not so far in G')
                    # print('shortest',len(nodes_to_add))
                    edges_to_add = [(nodes_to_add[i], nodes_to_add[i + 1]) for i in range(len(nodes_to_add) - 1)]
                    # print(edges_to_add)
                    # print(source,target,len_)
                    edges = edges + edges_to_add

    return edges


def edge_correction_v2(G_bfs, G_pre_extracted, gbfs_threshold, gpe_threshold):
    deg = nx.degree_centrality(G_bfs)
    N = len(G_bfs.nodes)
    for node in deg.keys():
        deg[node] = round(deg[node] * (N - 1))

    print(N)
    leaves = [node for node in G_bfs.nodes if deg[node] == 1]
    print(len(leaves))
    edges = []
    source_target_pair = []

    # leaves = leaves[200:300]

    for i in range(len(leaves) - 2):
        print('i', i)
        source = leaves[i]
        for j in range(i, len(leaves) - 1):
            target = leaves[j]

            len_ = nx.shortest_path_length(G_bfs, source=source, target=target)

            if len_ > gbfs_threshold:
                # print('=',source,target,'l',len_)
                # print('1.something to add!')
                nodes_to_add = nx.shortest_path(G_pre_extracted,
                                                source=source,
                                                target=target)
                if len(nodes_to_add) < gpe_threshold:
                    # N+=1
                    # print('2.not so far in G')
                    # print('shortest',len(nodes_to_add))
                    edges_to_add = [(nodes_to_add[i], nodes_to_add[i + 1]) for i in range(len(nodes_to_add) - 1)]
                    # print(edges_to_add)
                    print(source, target, len_)
                    edges = edges + edges_to_add
                    source_target_pair.append([source, target])

    return edges, source_target_pair


def filtering(Graph,
              beta_d,
              min_,
            folder_name,
            terminal_info,
            weighting_method_simplification='ER',
            BP_weights='BPtdens',
            weight_flag='length',
              input_flag = None):
    '''

    :param Graph: a networkx graph to be filtered.
    :param beta_d: beta input for the DMK solver.
    :param min_: threshold for the weights of the edges after filtering.
    :param folder_name: folder path where the outputs will be stored. It should be written "./runs/folder_name".
    :param terminal_info: for dat files (i.e. from continuous): [
            source_flag,
              sink_flag,
              btns_factor_source,
                btns_factor_sink].
                for images:
    terminal_info = [
            terminals,
            entries]
    :param weighting_method_simplification: 'ER', 'IBP', 'BPW'.
    :param BP_weights: 'BPtdens' to use optimal transport density as weights for edges, 'BPflux' to use optimal flux.
    :param weight_flag: 'length', to use the length of the edges; else, to use unit length edges.
    :param input_flag: 'image' or None (for dat files)
    :return:
        G_filtered: filtered graph (networkx graph).
        newGraph: dictionary, s.t., newGraph[i]= i-th cc (labeled from 0 to len(cc)-1).
        ncc: number of connected components of Graph.
        possible_terminals_source: for each i, possible_terminals_source[i]= "sources" of i-th cc.
        possible_terminals_sink: for each i, possible_terminals_sink[i]= "sources" of i-th cc.
        mapping: for each i, mapping[i]: labels of i-th cc -------------> labels of Graph.
    '''
    # ------------------------------------------------ filtering -------------------------------------------
    # Init dicts
    if len(terminal_info) == 2:# this is for image processing
        terminals = terminal_info[0]
        entries = terminal_info[1]
    elif len(terminal_info)==4:
        source_flag= terminal_info[0]
        sink_flag = terminal_info[1]
        btns_factor_source= terminal_info[2]
        btns_factor_sink= terminal_info[3]
    else:
        print('invalid terminal_info input!')

    conv_report = {}
    newGraph = {}
    mapping = {}
    inv_mapping = {}
    possible_terminals_source = {}
    possible_terminals_sink = {}
    edge_mapping = {}
    G_simplification = {}

    # Copy of Graph

    G = Graph.copy()

    #defining the beta_d for the simulations

    utils.updating_beta_discrete(beta_d)

    # Generating the cc-based graphs and the corresponding mappings
    newGraphList, mappingList, inv_mappingList,_ = utils.pickle2pygraph(G)

    # Iterating over the subgraphs
    ncc = len(newGraphList)
    for i in range(1, ncc + 1):

        conv_report[i]=[]
        newGraph[i] = newGraphList[i - 1]
        mapping[i] = mappingList[i - 1]
        inv_mapping[i] = inv_mappingList[i - 1]
        # ------------------------------------  filtering: this is done in a different way for images -----------------
        if len(terminal_info) == 2:
            source_candidates = [terminals[entry] for entry in entries]
            possible_terminals_source[i] = [inv_mapping[i][node] for node in source_candidates]
            possible_terminals_sink[i] = [inv_mapping[i][node] for node in terminals if node not in source_candidates]

        else:

            possible_terminals_source[i], possible_terminals_sink[i] = terminals_from_cont(newGraph[i], source_flag,
                                                                                       sink_flag,
                                                                                       btns_factor_source,
                                                                                       btns_factor_sink)

        edge_mapping[i] = utils.pygraph2dat(newGraph[i], possible_terminals_source[i], possible_terminals_sink[i], i,
                                      folder_name, mapping[i], input_flag)


        #print('executing graph2incidence_matrix for the component %s' % i)
        utils.using_graph2incidence_matrix(folder_name, i, weight_flag)
        print('_____________________________EXECUTING BP solver___________________________________________')
        BP_solver(folder_name, i)
        data_folder_name = folder_name + '/component' + str(i)
        G_simplification[i] = utils.dat2pygraph(newGraph[i], data_folder_name, edge_mapping[i],
                                              min_, BP_weights)

        # Checking convergence

        conver_bool = utils.convergence_tester( "../otp_utilities/muffe_sparse_optimization/simplifications/"+folder_name + '/component' + str(i))

        conv_report[i].append('|last_var_tdens|>1E-20:'+str(conver_bool))

        # Defining terminal labels for sources and sinks

        for node in G_simplification[i].nodes():
            node_in_original, _,_ = quality_measure.old_label(node, newGraph[i], G_simplification[i])
            if node_in_original in possible_terminals_source[i]:
                ss_label = 1
            elif node_in_original in possible_terminals_sink[i]:
                ss_label = -1
            else:
                ss_label = 0
            G_simplification[i].nodes[node]['terminal'] = ss_label

        isolated_nodes = len([node for node in G_simplification[i].nodes() if G_simplification[i].degree(node)==0])

        node_test = len(newGraph[i].nodes) - (len(G_simplification[i].nodes) - isolated_nodes)
        conv_report[i].append('|V(G)-V(G_f)|= '+ str(node_test))

        edge_test = len(newGraph[i].edges) - len(G_simplification[i].edges)
        conv_report[i].append(('|E(G)-E(G_f)|= '+ str(edge_test)))

        ccomp_test = len(list(nx.connected_components(G_simplification[i])))-isolated_nodes
        conv_report[i].append('|cc(G_f)| = '+ str(ccomp_test))

        test1 = conver_bool == True
        test2 =  node_test > 0
        test3 = edge_test > 0
        test4 = ccomp_test == 1 #fix this

        test_val = [test1,test2,test3,test4]

        conv_report[i]= ['success rate (%)= '+ str(int((sum(test_val)/4)*100))] + conv_report[i]+test_val
    # Assigning opt_pot to the nodes

    # Joining all the simplifications into a single graph

    #G = newGraphList[0]
    G_filtered = G_simplification[1]
    for i in range(1, ncc):
        #G = nx.disjoint_union(G, newGraphList[i])
        G_filtered = nx.disjoint_union(G_filtered, G_simplification[i + 1])

    # Adding the terminals (to track disconnectivities) (to do this we need to comment the nex line. But it's not working yet)

    G_filtered.remove_nodes_from(
        [x for x in G_filtered.nodes() if G_filtered.degree(x) == 0])


    # Reweighting the edges

    deg = nx.degree_centrality(G_filtered)
    posGsimpl = nx.get_node_attributes(G_filtered, 'pos')

    if weighting_method_simplification == 'ER':
        N = len(G_filtered.nodes())
        for edge in G_filtered.edges():
            if deg[edge[0]] * deg[edge[1]] != 0:
                G_filtered.edges[(edge[0], edge[1])]['tdens'] = G_filtered.nodes[edge[0]][
                                                                                 'weight'] / (
                                                                                     deg[edge[0]] * (N - 1)) + \
                                                                             G_filtered.nodes[edge[1]][
                                                                                 'weight'] / (
                                                                                     deg[edge[1]] * (N - 1))
                

    elif weighting_method_simplification == 'IBP':
        print('relabel/reweig!')
        G_filtered_relabeled = relabeling(G_filtered, G)
        G_filtered = reweighting(G_filtered_relabeled, G)
        # posGsimpl = nx.get_node_attributes(G_filtered, 'pos')
    elif weighting_method_simplification == 'BPW':
        print('BPW used!')
        pass

    print('number of edges after reweighting:',len(G_filtered.edges()))

    return G_filtered, newGraph, ncc, possible_terminals_source, possible_terminals_sink, mapping, conv_report

    # ------------ end of filtering -----------------------------------------


def img2filtering(image_path, new_size, number_of_colors, t1, t2, number_of_cc, graph_type, beta_d,weighting_method_simplification='ER',entries=[0]):
    '''
    This takes as input an image and outputs the filtered graph.
    :param image_path: string.
    :param number_of_colors: number of colors for the output image.
    :param t1: noise threshold. If t=0, then all the new pixels are preserved with their new colors.
    :param t2: threshold for the weights of the edges after pre-extraction.
    :param number_of_cc: number of connected components of the graph represented by the image. If None, then only 1
    cc is assumed.
    :param graph_type: 1 (to use edges and vertices of the grid), 2 (to use only edges).
    :param beta_d: beta input for the DMK solver.
    :return:
        G_filtered: filtered graph (networkx graph).
    '''

    # reading the image, doing pre extraction, getting bfs approx
    pre_extraction.bfs_preprocess(image_path,new_size, number_of_colors, t1, t2, number_of_cc, graph_type)
    #filtering the bfs graph
    G_filtered = img_pre_extr2filtering(image_path, filter_size, beta_d,weighting_method_simplification,entries)


    return G_filtered


#-------------------------test1------------------------------------------------------

filter_size=0.045
image_path = "./runs/graph_from_image/image.jpg"
weighting_method_simplification = 'ER'
beta_d=1.0
new_size = 100

#img_pre_extr2filtering(image_path, filter_size, beta_d,weighting_method_simplification,entries)

#--------------------------test2----------------------------------------------------

new_size=100
ratio=new_size/1200
#print('ratio:',ratio)
t1=0.1
t2=.5
image_path = "./runs/graph_from_image/image.jpg"
number_of_cc=1
number_of_colors=50
graph_type='1'

#img2filtering(image_path, new_size, number_of_colors, t1, t2, number_of_cc, graph_type, beta_d)
