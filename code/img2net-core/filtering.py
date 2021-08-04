## import stuff
import networkx as nx
import numpy as np
import sys
import pickle as pkl

from scipy.spatial import ConvexHull, convex_hull_plot_2d
from scipy.spatial import distance

# Accesing root path
import os

file_path = os.path.dirname(os.path.realpath(__file__))
with open("/Users/dtheuerkauf/Nextrout/nextrout_location.txt") as f:
    lines = f.readlines()
root = lines[0]
print('root',root)

# Import I/O for timedata
try:
    sys.path.append(root + "/../dmk/globals/python/timedata/")
    import timedata as td
except:
    print("Global repo non found")

# Import geometry tools
sys.path.append(root + "/../dmk/geometry/python/")
import meshtools as mt

sys.path.append(root + "/../dmk/dmk_solver/otp_solver/preprocess/assembly/")
import example_grid

# Import dmk tools
sys.path.append(root + "/../dmk/dmk_solver/otp_solver/python/")
import dmk_p1p0

sys.path.append(
    root + "/../dmk/dmk_solver/build/python/fortran_python_interface/"
)
from dmk import (
    Dmkcontrols,  # controls for dmk simulations)
    Timefunctionals,  # information of time/algorithm evolution
    Dmkinputsdata,  # structure variable containg inputs data
    build_subgrid_rhs,  # procedure to preprocess forcing term f
    Tdenspotentialsystem,  # input/output result tdens, pot
    dmkp1p0_steady_data,  # main interface subroutine
)

# Import plot tools
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

sys.path.append(root + "/../dmk/dmk_solver/graph_otp_solver/python")
import dmk_graph


def concatenate(lists):
    """
    This concatenates all the lists contained in a list.
    :param lists: a list of lists.
    :return:
        new_list: concatenated list.
    """
    new_list = []
    for i in lists:
        new_list.extend(i)
    return new_list


def terminals_from_cont(
    Graph,
    forcing_flag,
    extra_info,
    btns_factor_source,
    btns_factor_sink,
    terminal_criterion="branch_convex_hull+btns_centr",
):
    """
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
    """

    # Defining source-sink nodes of the graph (nodes inside the source or sink of the continuous DMK)

    nodes_in_source = []
    nodes_in_sink = []

    for node in Graph.nodes():
        terminal_val = Graph.nodes[node]["terminal"]
        if terminal_val == 1:
            nodes_in_source.append(node)
        elif terminal_val == -1:
            nodes_in_sink.append(node)

    
    bn = nx.betweenness_centrality(Graph, normalized=True)

    # min bn inside the source and sink

    max_bn_source = max([bn[node] for node in nodes_in_source])
    max_bn_sink = max([bn[node] for node in nodes_in_sink])

    # Defining source-sink candidates based only on btns

    kind_of_leaf_nodes_source = [
        key for key in nodes_in_source if bn[key] <= btns_factor_source * max_bn_source
    ]  # *(min_bn_source)+.0001]
    kind_of_leaf_nodes_sink = [
        key for key in nodes_in_sink if bn[key] <= btns_factor_sink * max_bn_sink
    ]  # *(min_bn_sink+.0001)]

    # Defining the subgraphs induced by the candidates

    sub_source = Graph.subgraph(kind_of_leaf_nodes_source)
    sub_sink = Graph.subgraph(kind_of_leaf_nodes_sink)

    # Removing repeated labels

    possible_terminals_source = set(kind_of_leaf_nodes_source)
    possible_terminals_sink = set(kind_of_leaf_nodes_sink)

    # print('poss',possible_terminals_source,possible_terminals_sink)

    if terminal_criterion == "single":

        terminals = list(possible_terminals_source) + list(possible_terminals_sink)

        # possible_terminals_source = list(possible_terminals_source)

        extra_nodes = []

        for i in range(len(terminals) - 1):  # range(len(possible_terminals_source)-1):

            node1 = terminals[i]

            pos1 = Graph.nodes[node1]["pos"]

            stop = False

            for j in range(i + 1, len(terminals)):

                node2 = terminals[j]

                if node1 != node2 and not stop:

                    pos2 = Graph.nodes[node2]["pos"]
                    dst = distance.euclidean(pos1, pos2)

                    if dst < 0.1:

                        # check distance

                        if node1 not in extra_nodes:

                            stop = True

                            extra_nodes.append(node1)

        possible_terminals_source = [
            int(node) for node in possible_terminals_source if node not in extra_nodes
        ]
        possible_terminals_sink = [
            int(node) for node in possible_terminals_sink if node not in extra_nodes
        ]

    elif terminal_criterion != "btns_centr":
        # Getting the coordinates to compute convex hulls

        coordinates_in_source = np.asarray(
            [
                [Graph.nodes[node]["pos"][0], Graph.nodes[node]["pos"][1]]
                for node in nodes_in_source
            ]
        )
        coordinates_in_source_list = concatenate(list(coordinates_in_source))

        coordinates_in_sink = np.asarray(
            [
                [Graph.nodes[node]["pos"][0], Graph.nodes[node]["pos"][1]]
                for node in nodes_in_sink
            ]
        )
        coordinates_in_sink_list = concatenate(list(coordinates_in_sink))

        # If the number of coordinates (or nodes) is not more than 3, then no convex hull computation

        if len(coordinates_in_source) >= 4 and len(coordinates_in_sink) >= 4:

            # Computing convex hull for the branches

            if terminal_criterion == "branch_convex_hull+btns_centr":
                source_hull = ConvexHull(coordinates_in_source)
                index_source_hull = np.asarray(source_hull.vertices)
                nodes_source_hull = []
                coord_source_hull = np.asarray(
                    [coordinates_in_source[node] for node in index_source_hull]
                )
                for node in nodes_in_source:
                    x, y = Graph.nodes[node]["pos"]
                    if x in coord_source_hull[:, 0] and y in coord_source_hull[:, 1]:
                        nodes_source_hull.append(node)

                sink_hull = ConvexHull(coordinates_in_sink)
                index_sink_hull = np.asarray(sink_hull.vertices)
                nodes_sink_hull = []
                coord_sink_hull = np.asarray(
                    [coordinates_in_sink[node] for node in index_sink_hull]
                )
                for node in nodes_in_sink:
                    x, y = Graph.nodes[node]["pos"]
                    if x in coord_sink_hull[:, 0] and y in coord_sink_hull[:, 1]:
                        nodes_sink_hull.append(node)
                possible_terminals_source = set(
                    kind_of_leaf_nodes_source + nodes_source_hull
                )
                possible_terminals_sink = set(kind_of_leaf_nodes_sink + nodes_sink_hull)

            # Computing convex hull for all the nodes defined as candidates

            elif terminal_criterion == "whole_convex_hull+btns_centr":  # not working!
                single_source_hull = ConvexHull(coordinates_in_source_list)
                single_index_source_hull = np.asarray(single_source_hull.vertices)
                nodes_source_hull = []
                single_coord_source_hull = np.asarray(
                    [
                        coordinates_in_source_list[node]
                        for node in single_index_source_hull
                    ]
                )
                for node in nodes_in_source:
                    x, y = Graph.nodes[node]["pos"]
                    if (
                        x in single_coord_source_hull[:, 0]
                        and y in single_coord_source_hull[:, 1]
                    ):
                        nodes_source_hull.append(node)
                single_sink_hull = ConvexHull(coordinates_in_sink_list)
                single_index_sink_hull = np.asarray(single_sink_hull.vertices)
                nodes_sink_hull = []
                single_coord_sink_hull = np.asarray(
                    [coordinates_in_sink_list[node] for node in single_index_sink_hull]
                )
                for node in nodes_in_sink:
                    x, y = Graph.nodes[node]["pos"]
                    if (
                        x in single_coord_sink_hull[:, 0]
                        and y in single_coord_sink_hull[:, 1]
                    ):
                        nodes_sink_hull.append(node)
                possible_terminals_source = set(
                    kind_of_leaf_nodes_source + nodes_source_hull
                )
                possible_terminals_sink = set(kind_of_leaf_nodes_sink + nodes_sink_hull)

    return possible_terminals_source, possible_terminals_sink


def filtering(
    Gpe,
    sources=None,
    sinks=None,
    beta_d=1.5,
    threshold=1e-3,
    tdens0=None,
    BPweights="tdens",
    stopping_threshold_f=1e-3,
    weight_flag="unit",
    rhs=None,
    MaxNumIter = 100,
    verbose=False,
):

    inputs = {}

    if sources is None and sinks is None and rhs is None:

        raise ValueError("Either rhs or sources/sinks need to be passed as inputs.")

    ### relabeling

    # todo: add an if for the case in which nodes are already relabeled

    mapping = {}
    k = -1
    for node in Gpe.nodes():
        k += 1
        mapping[node] = k
    Gpe_rel = nx.relabel_nodes(Gpe, mapping, copy=True)

    edges = Gpe_rel.edges()
    nedges = len(edges)
    nodes = Gpe_rel.nodes()
    nnodes = len(nodes)

    # tdens0

    if tdens0 != None:
        try:
            tdens0 = np.array([(Gpe_rel.edges[edge]["tdens"]) for edge in edges])
        except:
            tdens0 = np.array([(Gpe_rel.edges[edge]["flux"]) for edge in edges])

    # topol

    topol = np.zeros((nedges, 2))
    k = -1
    for edge in edges:
        k += 1
        topol[k, :] = edge

    # weight (uniform)

    weight = np.empty(nedges, dtype=object)

    k = -1
    for edge in edges:
        k += 1
        if weight_flag == "unit":
            weight[k] = 1
        elif weight_flag == "length":
            weight[k] = distance.euclidean(
                Gpe_rel.nodes[edge[0]]["pos"], Gpe_rel.nodes[edge[1]]["pos"]
            )
        else:
            weight[k] = Gpe_rel.edges[edge][weight_flag]

    # rhs (f+ and f-)

    if (
        sinks is not None and sources is not None
    ):  # there are lists from the sources and sinks are going to be chosen.
        # (else) if this is not pass, then the rhs is passed.

        rhs = np.zeros(nnodes)
        sources_rel = [mapping[node] for node in sources]
        sinks_rel = [mapping[node] for node in sinks]

        number_sources = len(sources_rel)
        number_sinks = len(sinks_rel)

        for node in nodes:
            if node in sources_rel:
                rhs[node] = 1 / number_sources
            elif node in sinks_rel:
                rhs[node] = -1 / number_sinks
            else:
                rhs[node] = 0
    else:
        sources_rel = [i for i in range(len(rhs)) if rhs[i] > 0]
        sinks_rel = [i for i in range(len(rhs)) if rhs[i] < 0]

    
    assert sum(rhs) < 0.01
    assert len(rhs) == nnodes
    # init and set controls
    ctrl = Dmkcontrols.DmkCtrl()
    Dmkcontrols.get_from_file(ctrl, root + "/nextrout_core/dmk_discr.ctrl")
    # if and where save data
    ctrl.id_save_dat = 1
    ctrl.fn_tdens = "tdens.dat"
    ctrl.fn_pot = "pot.dat"
    ctrl.max_time_iterations = MaxNumIter
    # if and where save log
    ctrl.id_save_statistics = 1
    ctrl.fn_statistics = "dmk.log"
    # if print info
    #
    if verbose:
        ctrl.info_state = 3
        ctrl.info_update = 3
        print(ctrl.outer_solver_approach)
    else:
        ctrl.info_state = 0
        ctrl.info_update = 0

    [info, tdens, pot, flux, timefun] = dmk_graph.dmk_graph(
        topol,
        rhs,
        pflux=beta_d,
        tdens0=tdens0,
        tolerance = stopping_threshold_f,
        weight=weight,
        ctrl=ctrl,
    )

    tdens = list(tdens)
    flux = list(flux)

    if (info == 0) and verbose:
        print("Convergence achieved")

    max_flux = max(flux)
    max_tdens = max(tdens)
    Gf = nx.Graph()
    ed_count = -1
    weights_in_Gf = []
    for edge in Gpe_rel.edges():
        ed_count += 1
        if BPweights == "flux":
            if abs(flux[ed_count]) > max_flux * threshold:
                Gf.add_edge(*edge, flux=flux[ed_count])
                weights_in_Gf.append(flux[ed_count])

        elif BPweights == "tdens":
            if abs(tdens[ed_count]) > max_tdens * threshold:
                Gf.add_edge(*edge, tdens=tdens[ed_count])
                weights_in_Gf.append(tdens[ed_count])

        else:
            raise ValueError("BPweights flag not defined!.")
        try:
            Gf.add_node(
                edge[0], weight=Gpe_rel.nodes[edge[0]]["tdens"]
            )  # todo: this needs to be fixed once the flux is working again (BPweights)
            Gf.add_node(edge[1], weight=Gpe_rel.nodes[edge[1]]["tdens"])
        except:
            pass

    Gf.remove_nodes_from(list(nx.isolates(Gf)))

    weights_in_Gf = np.array(weights_in_Gf)
    colors = []

    for node in Gf.nodes():

        Gf.nodes[node]["pos"] = Gpe_rel.nodes[node]["pos"]

        if node in sources_rel:
            colors.append("g")
        elif node in sinks_rel:
            colors.append("r")
        else:
            colors.append("k")

    inputs["topol"] = topol
    inputs["rhs"] = rhs
    inputs["pflux"] = beta_d
    inputs["tdens0"] = tdens0

    return Gf, weights_in_Gf, colors, inputs


def bifurcation_paths(G, terminals):
    """
    This script takes a filtered graph and reduces its paths (sequences of nodes with degree 2) to a single edge.

    :param G:  filtered graph (networkx graph).
    :param terminals: union of source and sink nodes.
    :return:
        G: reduced graph.
    """

    G = G.copy()
    N = len(G.nodes())
    deg_norm = nx.degree_centrality(G)
    deg = {}
    for key in deg_norm.keys():
        deg[key] = int(round((N - 1) * deg_norm[key]))
    # print(deg)
    # deg_neq_2=[node for node in G.nodes() if deg[node]!= 2]
    deg_3 = [
        node
        for node in G.nodes()
        if deg[node] >= 3 or deg[node] == 1 or node in terminals
    ]

    G_wo_bifuc = G.copy()
    for node in deg_3:
        G_wo_bifuc.remove_node(node)
    cc = list(nx.connected_components(G_wo_bifuc))
    # print(list(cc))
    connect_points = {}
    index = 1
    for comp in cc:
        comp = set(comp)
        neighs = {
            neigh for node in comp for neigh in G.neighbors(node) if neigh not in comp
        }
        # print(neighs)
        assert len(neighs) == 2
        G.remove_nodes_from(comp)
        G.add_edge(*tuple(neighs))

    return G




def filtering_from_image(small_G_filtered, sources, sinks, beta_d,  color_dict, partition_dict, weighting_method_simplification='ER',folder_name='./runs/no_name/',plotting=True):
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

    G_filtered,_,_,_ = filtering(
                            Graph,
                            sources=sources,
                            sinks=sinks,
                            beta_d = beta_d,
                            threshold = min_,
                            BPweights="tdens",
                            stopping_threshold_f=1e-5,
                            weight_flag="length"
                        )

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

        for node in small_G_filtered.nodes():
            if node in sources:
                color = 'green'
            elif node in sinks:
                color = 'red'
            else:
                color = 'black'
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

        for node in small_G_filtered.nodes():
            if node in sources:
                color = 'green'
            elif node in sinks:
                color = 'red'
            else:
                color = 'black'
            x = small_G_filtered.nodes[node]['pos'][0]
            y = small_G_filtered.nodes[node]['pos'][1]
            circle1 = plt.Circle((x, y), .015, color=color, fill=False, lw=4)
            ax.add_artist(circle1)

        plt.savefig(folder_name + '/filtered_graph+terminals', transparent=False)

    with open(folder_name + '/filtered_graph.pkl', 'wb') as file:
        pkl.dump(G_filtered, file)
    plt.close('all')

    return G_filtered



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
