import networkx as nx
import numpy as np
import itertools
import shutil
import matplotlib.pyplot as plt
import pickle as pkl
import os
import time
import decimal
import cv2 as cv
from PIL import Image
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import time
#--------------------------------------------
import quality_measure
#--------------------------------------------

def prod_dict(File1, File2, weights1, weights2):
	'''
	Key-wise multiplication of two dictionaries. Component-wise multiplication of two lists.
	:param File1:  dictionary.
	:param File2: dictionary.
	:param weights1: list.
	:param weights2: list.
	:return:
		newDict: dictionary with same keys as File1 (or File2) but values equal to the multiplication of File1 and File2.
		newWeights: list whose entries are the multiplication of weights1 and weights2.
	'''

	newDict = {}
	newWeights = []
	if File1.keys() != File2.keys():
		print("not the same keys!")
	else:
		i = 0
		for key in File1.keys():
			newDict[key] = File1[key] * File2[key]
			newWeights.append(weights1[i] * weights2[i])
			i += 1
		return newDict, newWeights


def dict2graph(bar_pos, dict_weights_func):
	'''
	Graph generator from two dictionaries.
	:param bar_pos: dictionary, s.t., dictionary[key]= position of key in domain.
	:param dict_weights_func: dictionary, s.t., dictionary[key]= weight of key.
	:return:
		X_func: weighted graph whose nodes have two attributes: 'pos' and 'weight'.
		pos: nx.get_node_attributes(X_func, "pos").
	'''
	X_func = nx.Graph()
	X_func.add_nodes_from(bar_pos.keys())
	# print(X_func.nodes())
	for n, p in bar_pos.items():
		X_func.nodes[n]["pos"] = p
		X_func.nodes[n]["weight"] = dict_weights_func[int(n)]
	pos = nx.get_node_attributes(X_func, "pos")

	return X_func, pos


def node_size(dict_weights_func):
	'''
	Defining the size of the nodes. This gives the size to the nodes depending on the weight of each one of them.
	Plot related.
	:param dict_weights_func: dictionary, s.t., dictionary[key]= weight of key.
	:return:
		size_func: list whose entries are the sizes for the nodes.
	'''
	size_func = [
		dict_weights_func[v] * 10 for v in dict(sorted(dict_weights_func.items()))
	]
	return size_func

# not sure if used
def plot_graph(X_func, size_func, file_name):
	pos_func = nx.get_node_attributes(X_func, "pos")
	if file_name == "source.dat":
		color = "b"
	elif file_name == "sink.dat":
		color = "g"
	else:
		color = "r"
	nx.draw_networkx(X_func, pos=pos_func, node_color=color, node_size=size_func)
	# print(folder_name+'/graph_pc'+str(int(min_*100))+'_graph'+graph_type+'.dat')
	# plt.savefig(folder_name+'/graph_pc'+str(int(min_*100))+'_graph'+graph_type+'.png')
	plt.axis("on")
	plt.show(block=False)




def get_first_neig(node, dict_seq):
	'''
	This returns the vertices that belong to the triangle for which 'node' is the barycenter
	:param node: node in G_bar.
	:param dict_seq:  dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular and squared grids.
	:return:
		same_triang: all the nodes in the triangle (or element of the grid) s.t., node is its barycenter.
	'''
	same_triang = list(dict_seq[node])
	return same_triang



def get_sec_neig(node, dict_seq):
	'''
	This returns all the triangles that share either a vertex or an edge with the triangle in which the 'node'
	is the barycenter.
	:param node: target node.
	:param dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular
		and squared grids.
	:return:
		dict_sec_neig: dictionary, s.t., dict_sec_neig[node]= indexes of all the surrounding triangles of 'node'.
		index_: **check this output. it seems to be removable.**
	'''
	same_triang = get_first_neig(
		node, dict_seq
	)  # getting the nodes of the triang that has 'node' as barycenter
	dict_sec_neig = (
		{}
	)  # dict that contains all the surrounding triangles for the 'node'
	dict_sec_neig[node] = []

	# indexes of all the surrounding triangles. This will be useful to call not only the mentioned triangles but also their bar
	index_ = {}
	index_[node] = []
	for key in dict_seq.keys():
		# if one of the nodes is in the triangle, then ...
		if (
				same_triang[0] in dict_seq[key]
				or same_triang[1] in dict_seq[key]
				or same_triang[2] in dict_seq[key]
		):

			# (to avoid repeating the triangles)
			repeated_tr = list(np.array(same_triang) == np.array(dict_seq[key]))

			if repeated_tr != [True, True, True]:
				# ... we define the triangle indexes for a particular node
				dict_sec_neig[node].append(dict_seq[key])
				index_[node].append(key)
	return dict_sec_neig, index_



def get_sec_neig_edges(node, dict_seq):
	'''
	This returns all the triangles that share an edge with the triangle in which the 'node' is the barycenter
	:param node: target node.
	:param dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular
		and squared grids.
	:return:
		dict_sec_neig: dictionary, s.t., dict_sec_neig[node]= indexes of all the surrounding triangles of 'node'.
		index_: **check this output. it seems to be removable.**
	'''
	same_triang = get_first_neig(
		node, dict_seq
	)  # getting the nodes of the triangle that has 'node' as barycenter
	dict_sec_neig = (
		{}
	)  # dict that contains all the surrounding triangles for the 'node'
	dict_sec_neig[node] = []
	index_ = (
		{}
	)  # indexes of all the surrounding triangles. This will be useful to call not only the mentioned triangles
	# but also their bar
	index_[node] = []
	for key in dict_seq.keys():
		if (
				(same_triang[0] in dict_seq[key] and same_triang[1] in dict_seq[key])
				or (same_triang[1] in dict_seq[key] and same_triang[2] in dict_seq[key])
				or (same_triang[2] in dict_seq[key] and same_triang[0] in dict_seq[key])
		):
			repeated_tr = list(np.array(same_triang) == np.array(dict_seq[key]))
			if repeated_tr != [True, True, True]:
				dict_sec_neig[node].append(dict_seq[key])
				index_[node].append(key)
	return dict_sec_neig, index_

def get_sec_neig_edges_square(node, dict_seq, node2box_index):
	'''
	This returns all the triangles that share an edge with the triangle in which the 'node' is the barycenter
	:param node: target node.
	:param dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular
		and squared grids.
	:return:
		dict_sec_neig: dictionary, s.t., dict_sec_neig[node]= indexes of all the surrounding triangles of 'node'.
		index_: **check this output. it seems to be removable.**
	'''
	same_triang = get_first_neig(node, dict_seq)  # getting the nodes of the triang that has 'node' as barycenter
	neighboring_boxes = []
	pair_of_vertices_same_triang = list(itertools.combinations(same_triang, 2))
	for pair in pair_of_vertices_same_triang:
		node1=pair[0]
		node2=pair[1]
		indexes_node1 = node2box_index[node1]
		indexes_node2 = node2box_index[node2]
		print(indexes_node1,indexes_node2)
		intersect = [index for index in indexes_node1 if index in indexes_node2]
		neighboring_boxes = neighboring_boxes + intersect
	neighboring_boxes = list(set(neighboring_boxes))
	neighboring_boxes.remove(node)
	# print(neighboring_boxes)
	dict_sec_neig = {}  # dict that contains all the surrounding triangles for the 'node'
	dict_sec_neig[node] = []

	# indexes of all the surrounding triangles. This will be useful to call not only the mentioned triangles but also their bar
	index_ = neighboring_boxes

	return dict_sec_neig, index_







def connecting_edges(G_bar, node, min_, graph_type, dict_seq, max_, weighting_method, input_flag=None,node2box_index=None, t3 = 1):
	'''
	Testing the condition of the tdens for a single barycenter . If condition is satisfied, add the edge with weight
	equal to the min. In-place modification of G_bar.
	:param G_bar: a weighted networkX graph whose nodes are the barycenters of the elements of the grid. No edges.
	:param node: target node.
	:param min_: threshold for the weights of the edges after pre-extraction.
	:param graph_type: 1 (to use edges and vertices of the grid), 2 (to use only edges), 3 (to use elements of the grid).
	:param dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular and squared grids.
	:param max_: maximum weight of nodes in G_bar.
	:param weighting_method: 'ER', 'AVG'.
	:param input_flag: 'image' or None (for dat files)
	:param node2box_index:  given a node, node2box_index[node] returns the indices of all the elements of the partition s.t.,
	node is a vertex of these elements
	:return:
	'''
	if graph_type == "1":
		# If True, then just the 'first neighbors' are taken.
		if input_flag != 'image':
			node = str(node)
			dict_sec_neig, index_ = get_sec_neig(node, dict_seq)
			index_ = index_[node]
		else:
			node=node-1
			dict_sec_neig, index_ = get_sec_neig_square(node, dict_seq, node2box_index)

	elif graph_type == "2":
		# If True, then just the 'second neighbors' are taken (firsts included).
		if input_flag != 'image':
			node = str(node)
			dict_sec_neig, index_ = get_sec_neig_edges(node, dict_seq)
			index_ = index_[node]
		else:
			node = node - 1
			dict_sec_neig, index_ = get_sec_neig_edges_square(node, dict_seq, node2box_index)


	#print(index_)
	for bar_ in index_:
		if graph_type == "1" or graph_type == "2":  # both edges and vertices
			# Checking if the weights are greater than the threshold x max
			if (
					 G_bar.nodes[bar_]["weight"] > min_ * max_
					and G_bar.nodes[node]["weight"] > min_ * max_
			):
				'''
				max_ * t3> G_bar.nodes[bar_]["weight"] > min_ * max_
						and max_*t3 > G_bar.nodes[node]["weight"] > min_ * max_
				'''
				if weighting_method == "AVG":
					# If True, then we assign weights to the edges based on the 'AVG' method
					w_list = [G_bar.nodes[bar_]["weight"], G_bar.nodes[node]["weight"]]
					w = sum(w_list) / len(w_list)  # <--- the average
					# w=w0 #<--- this is the minimum
					G_bar.add_edge(node, bar_, weight=w)
				elif weighting_method == "ER":
					# If True, then we just add the edge. Later on we define the weight.
					G_bar.add_edge(node, bar_)


def grid_filter(G_bar, G_triang, min_, dict_seq):
	'''
	This script adds the edges of a triangle to the graph type 3 if the weight of its barycenter is greater than the
	threshold (min_) x max_.
	:param G_bar: a weighted networkX graph whose nodes are the barycenters of the elements of the grid. No edges.
	:param G_triang: the grid graph.
	:param min_: threshold for the weights of the edges after pre-extraction.
	:param dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular and squared grids.
	:return:
		 G_pre_extracted: the input graph but with edges generated according to graph_type 3.
	'''

	nodes_dict = G_bar.nodes(data='weight')
	max_ = max([entry[1] for entry in nodes_dict])

	G_pre_extracted = G_triang.copy()

	# Getting the max key (why?)

	#keys = [int(key) for key in G_pre_extracted.nodes().keys()]
	#max_label = max(keys)

	# Removing the edges of the triangulation graph. We are left only with the nodes. Also the useless ones.

	edges_triang = list(G_pre_extracted.edges())
	G_pre_extracted.remove_edges_from(edges_triang)

	for bar_ in G_bar.nodes():

		# Getting the edges of the triangle for a given bar 'bar_'

		same_triang = get_first_neig(bar_, dict_seq)
		w = G_bar.nodes[bar_]["weight"]

		if w >= min_ * max_:
			# If True, we add the edges to the graph
			for node in same_triang:
				G_pre_extracted.nodes[node]['weight']=0
			edges = list(itertools.combinations(same_triang, 2))
			membership = [
				edges[0] in G_pre_extracted.edges(),
				edges[1] in G_pre_extracted.edges(),
				edges[2] in G_pre_extracted.edges(),
			]

			# We iterate over them to define their weights
			for i in [0, 1, 2]:
				bool_val = membership[i]
				# If the edge is in the graph (thanks to another barycenter), we sum the weight of the barycenter 'bar_' to the existing weight
				if bool_val == True:
					new_weight = (
							G_pre_extracted.edges[(edges[i][0], edges[i][1])]["weight"]
							+ float(w) / 2.0
					)
					G_pre_extracted.edges[(edges[i][0], edges[i][1])]["weight"] = new_weight
				# If the edge is not in the graph, we add half of the weight of the barycenter
				else:
					new_weight = float(w) / 2.0
					G_pre_extracted.add_edge(edges[i][0], edges[i][1], weight=new_weight)

	return G_pre_extracted


def node_edge_filter(G_bar, min_, graph_type, dict_seq, weighting_method,input_flag=None,node2box_index=None, t3 = 1):
	'''
	:param G_bar: a weighted networkX graph whose nodes are the barycenters of the elements of the grid. No edges.
	:param min_:  threshold for the weights of the edges after pre-extraction.
	:param graph_type: 1 (to use edges and vertices of the grid), 2 (to use only edges).
	:param dict_seq:  dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular and squared grids.
	:param weighting_method: 'ER', 'AVG'.
	:param input_flag: 'image' or None (for dat files)
	:param node2box_index: given a node, node2box_index[node] returns the indices of all the elements of the partition s.t.,
	node is a vertex of these elements
	:return:
		G_bar: the input graph but with edges generated according to graph_type.
	'''
	"""
	This script generates graph type 1 and 2.
	"""
	nodes_dict = G_bar.nodes(data='weight')
	max_ = max([entry[1] for entry in nodes_dict])
	reduced_node_list = [int(node) for node in G_bar.nodes() if G_bar.nodes[node]['weight']>min_ * max_]

	print('Barycenters whose weight is over the threshold:',len(reduced_node_list))
	# Iterate over all the numbers (<--> nodes) to test the condition about the threshold:
	for n in reduced_node_list:#range(len(G_bar.nodes()))
		connecting_edges(
			G_bar, n + 1, min_, graph_type, dict_seq, max_, weighting_method,input_flag,node2box_index, t3 = t3
		)

	if weighting_method == "ER":

		# Compute degree centrality

		deg = nx.degree_centrality(G_bar)
		N = len(G_bar.nodes())

		# If True, then apply 'ER' method

		for edge in G_bar.edges():
			G_bar.edges[(edge[0], edge[1])]["weight"] = G_bar.nodes[edge[0]]["weight"] / (
					deg[edge[0]] * (N - 1)
			) + G_bar.nodes[edge[1]]["weight"] / (deg[edge[1]] * (N - 1))
	elif weighting_method == "AVG":
		pass
	else:
		print("weighting_method not defined")

	return G_bar


"""
Filtering triangles
"""


def pre_extraction(
		G_bar, G_triang, dict_seq, min_, graph_type='1',weighting_method='ER'
):
	'''
	:param G_bar: a weighted networkX graph whose nodes are the barycenters of the elements of the grid.
	:param G_triang: the grid graph.
	:param dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular and squared grids.
	:param graph_type: 1 (to use edges and vertices of the grid), 2 (to use only edges), 3 (to use elements of the grid).
	:param min_: threshold for the weights of the edges after pre-extraction.
	:param max_: maximum weight of nodes in G_bar.
	:param weighting_method: 'ER', 'AVG'.
	:return:
		G_pre_extracted: pre-extracted graph.
	'''


	if graph_type == "3" and weighting_method == "ER":
		print("this wm does not apply for graph type 3. Try with 'AVG'.")
	else:

		if graph_type == "1" or graph_type == "2":
			G_pre_extracted = G_bar.copy()
			# filtering
			edges_ = list(G_pre_extracted.edges())
			G_pre_extracted.remove_edges_from(edges_)
			G_pre_extracted = node_edge_filter(
				G_pre_extracted, min_, graph_type, dict_seq, weighting_method
			)
		elif graph_type == "3":
			# filtering
			# print('here')
			G_pre_extracted = grid_filter(G_bar, G_triang, min_, dict_seq)

		"""
		Removing the non-useful barycenters
		"""

		G_pre_extracted.remove_nodes_from(list(nx.isolates(G_pre_extracted)))

		return G_pre_extracted


def coloring(pixel_val, N):
	'''
	This functions assigns a "color" (an integer between 0 and N-1) to a pixel.
	:param pixel_val: normalized real value of the pixel.
	:param N: number of colors.
	:return:
		color: integer \in \{0,1,...,N-1\}.

	'''
	intervals = np.linspace(0, 1, N, endpoint=True)
	interval_bool = [intervals[i] <= pixel_val <= intervals[i + 1] for i in range(len(intervals) - 1)]
	color = interval_bool.index(True)
	return color




def get_sec_neig_square(node, dict_seq, node2box_index):
	'''
	This returns all the squares that share either a vertex or an edge with the square in which the 'node'
	is the barycenter.
	:param node: target node.
	:param dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular
		and squared grids.
	:return:
		dict_sec_neig: dictionary, s.t., dict_sec_neig[node]= indexes of all the surrounding triangles of 'node'.
		index_: **check this output. it seems to be removable.**
	'''
	same_triang = get_first_neig(node, dict_seq)  # getting the nodes of the triang that has 'node' as barycenter
	neighboring_boxes = []
	for vertex in same_triang:
		neighboring_boxes = neighboring_boxes + node2box_index[vertex]
	neighboring_boxes = list(set(neighboring_boxes))
	neighboring_boxes.remove(node)
	# print(neighboring_boxes)
	dict_sec_neig = {}  # dict that contains all the surrounding triangles for the 'node'
	dict_sec_neig[node] = []

	# indexes of all the surrounding triangles. This will be useful to call not only the mentioned triangles but also their bar
	index_ = neighboring_boxes

	return dict_sec_neig, index_

def bar_square(coord):
	'''
	This returns the coordinates of the barycenter of a square defined by coord.
	:param coord: list of coordinates of square (len(4)).
	:return:
		x_bar: x-coordinate of barycenter.
		y_bar: y-coordinate of barycenter.
	'''
	x1,y1 = coord[0]
	x2,y2 = coord[3]
	x_bar = .5*(x1+x2)
	y_bar = .5*(y1+y2)
	return x_bar,y_bar


def resizing_image(image_path, number_of_colors, new_size,t=0,reversed_colors = True, ds_interpolation = 'nearest', plotting = True):
	'''
	This resizes and repaints an image.
	:param image_path: string.
	:param number_of_colors: number of colors for the output image.
	:param t: noise threshold. If t=0, then all the new pixels are preserved with their new colors.
	:return:
		width: width of the output image.
		color_dict: dictionary s.t., color_dict[key]= real value for the key-th pixel.
		key is the index for the pixels in the resized image.
		saving_path: string, to where the new image is saved.
	'''
	if ds_interpolation == 'nearest':
		ds_funct = cv.INTER_NEAREST
	elif ds_interpolation == 'area':
		ds_funct = cv.INTER_AREA

	last_word = image_path.split('/')[-1]
	new_folder_name = last_word.split('.')[0]
	saving_path = './runs/'+new_folder_name

	try:
		os.mkdir(saving_path)
	except:
		pass
	#print('resizing',os.getcwd(), image_path)
	im = Image.open(image_path)
	width, height, = im.size
	pixel_values = np.array(list(im.getdata()))
	im = cv.imread(image_path)

	if new_size == 'automatic':
		new_size = max( int(0.30*width) , 100)

	if width != new_size:
		#resizing it

		img_rotate_90_clockwise = cv.rotate(im, cv.ROTATE_90_CLOCKWISE)

		img_ratio = new_size / width
		#print('ratio:', img_ratio)

		small_to_large_image_size_ratio = img_ratio
		small_img = cv.resize(img_rotate_90_clockwise,  # original image
							  (0, 0),  # set fx and fy, not the final size
							  fx=small_to_large_image_size_ratio,
							  fy=small_to_large_image_size_ratio,
							  interpolation=ds_funct)
		cv.imwrite(saving_path+'/resized_'+new_folder_name+'.jpg', small_img)
		res_im = Image.open(saving_path+'/resized_'+new_folder_name+'.jpg', 'r')
		width, height = res_im.size
		pixel_values = np.array(list(res_im.getdata()))
	else:
		print('new size = current size')


	partition_dict, dict_seq, node2box_index = quality_measure.partition_set(width + 1)

	number_of_colors += 1
	i = 0
	color_dict = {}
	for r, b, g in pixel_values:
		rgb = ((r & 0x0ff) << 16) | ((g & 0x0ff) << 8) | (b & 0x0ff)
		# print(rgb)
		color_dict[i] = int(rgb)
		i += 1
	#print(len(color_dict))
	colors = []
	max_ = max(color_dict.values())
	color_flag = 2

	for key in color_dict.keys():
		if reversed_colors == True:
			color_dict[key] = max(1 - color_dict[key] / max_, t)
		else:
			color_dict[key] = max(color_dict[key] / max_, t)
		# print(color_dict[key])
		#color = coloring(color_dict[key], number_of_colors)
		
		color = color_dict[key]

		colors.append(color)

	max_ = max(color_dict.values())


	if plotting:

		fig, ax = plt.subplots(1, 1, figsize=(17, 15))
		patches = []
		for key in partition_dict:
			square_edges = np.asarray([partition_dict[key][0]] + [partition_dict[key][2]] + [partition_dict[key][3]] + [
				partition_dict[key][1]] + [partition_dict[key][0]])

			s1 = Polygon(square_edges)
			patches.append(s1)
		p = PatchCollection(patches, alpha=1, linewidth=.0, edgecolor='b', cmap='Greys')

		p.set_array(np.array(colors))
		ax.add_collection(p)
		plt.colorbar(p)
		folder_path = image_path.replace(image_path.split('/')[-1],"")

		plt.savefig(saving_path+'/repainted_resized_image.png')
		#plt.show()

	return width,color_dict,saving_path

def weighted_partition2bar_graph(partition_dict, color_dict):
	'''
	This takes a partition and weights (colors) for its elements and outputs the G_bar.
	:param partition_dict:
	:param color_dict:
	:return:
		G_bar: a weighted networkX graph whose nodes are the barycenters of the elements of the grid. No edges.
	'''
	G_bar = nx.Graph()
	for key in partition_dict.keys():
		G_bar.add_node(key - 1)
		value = partition_dict[key]
		weight = color_dict[key - 1]
		G_bar.nodes[key - 1]['pos'] = np.array(bar_square(value))
		G_bar.nodes[key - 1]['weight'] = weight

	return G_bar



def pre_extraction_from_image(image_path,new_size,t2,number_of_colors=50,number_of_cc=1,graph_type='1',t1=0, reversed_colors = True, t3 = 1, ds_interpolation= 'nearest', plotting = False):
	'''
	This takes an image and return a graph extracted from it according to the pre-extraction rules.
	:param image_path: string.
	:param new_size: new size for the input image.
	:param t2: threshold for the weights of the edges after pre-extraction.
	:param number_of_colors: number of colors for the output image.
	:param number_of_cc: number of connected components of the graph represented by the image. If None, then only 1
	cc is assumed.
	:param graph_type: 1 (to use edges and vertices of the grid), 2 (to use only edges).
	:param t1: noise threshold. If t1=0, then all the new pixels are preserved with their new colors.
	:return:
		small_G_pre_extracted: pre-extracted graph.
	'''
	
	width, color_dict, folder_path = resizing_image(image_path, number_of_colors,new_size, t1, reversed_colors,ds_interpolation)

	N = width + 1
	_, G_triang = quality_measure.partition(N)
	partition_dict, dict_seq, node2box_index = quality_measure.partition_set(N)

	G_bar = weighted_partition2bar_graph(partition_dict, color_dict)

	#max_=max(color_dict.values())

	G_pre_extracted = G_bar.copy()

	# filtering
	edges_ = list(G_pre_extracted.edges())
	G_pre_extracted.remove_edges_from(edges_)

	G_pre_extracted = node_edge_filter(G_pre_extracted, t2, graph_type, dict_seq, 'ER','image',node2box_index, t3=t3)  # 12

	small_G_pre_extracted = G_pre_extracted.copy()
	small_G_pre_extracted.remove_nodes_from(list(nx.isolates(G_pre_extracted)))
	pos = nx.get_node_attributes(small_G_pre_extracted, 'pos')

	if plotting:

		fig, ax = plt.subplots(1, 1, figsize=(15, 15))
		patches = []
		for key in partition_dict:
			square_edges = np.asarray([partition_dict[key][0]] + [partition_dict[key][2]] + [partition_dict[key][3]] + [
				partition_dict[key][1]] + [partition_dict[key][0]])

			s1 = Polygon(square_edges)
			patches.append(s1)
		p = PatchCollection(patches, alpha=.7, cmap='YlOrRd', linewidth=.1, edgecolor='b')

		colors = np.array(list(color_dict.values()))
		p.set_array(colors)
		ax.add_collection(p)
		
		nx.draw_networkx(small_G_pre_extracted, pos, node_size=1, width=3, with_labels=False, edge_color='Gray', alpha=0.8,
						 node_color='black', ax=ax)

	if number_of_cc in [1,None] and plotting: #we assume then there's just one

		plt.savefig(folder_path + '/extracted_graph.png')

		cc_large = max(nx.connected_components(small_G_pre_extracted), key=len)
		small_G_pre_extracted = small_G_pre_extracted.subgraph(cc_large)

		pos = nx.get_node_attributes(small_G_pre_extracted, 'pos')
		nx.draw_networkx(small_G_pre_extracted, pos, node_size=3, width=3, with_labels=False, edge_color='black',
						 alpha=0.8, node_color='black', ax=ax)

	with open(folder_path + '/extracted_graph.pkl', 'wb') as file:
		pkl.dump(small_G_pre_extracted, file)

	with open(folder_path+'/real_image.pkl', 'wb') as file:
		pkl.dump(color_dict, file)

	with open(folder_path+'/G_bar.pkl', 'wb') as file:
		pkl.dump(G_bar, file)

	
	return small_G_pre_extracted, color_dict, partition_dict, folder_path


def tree_approximation(Graph, color_dict, partition_dict, folder_path, root = None):
	'''
	This returns a tree approximation of the input graph. In this case, the graph used is the bfs rooted at the lowest
	labeled node.
	:param Graph:  a networkx graph.
	:param root: node label to used as root for bfs.
	:return:
		bfs_Graph: bfs approximation of G.
	'''

	nodes = sorted(list(Graph.nodes))
	if root == None:
		root=nodes[0]
	else:
		pass

	bfs_Graph = nx.bfs_tree(Graph, root)
	bfs_Graph = bfs_Graph.to_undirected()

	for node in bfs_Graph.nodes():
		for word in ['pos', 'weight']:
			bfs_Graph.nodes[node][word] = Graph.nodes[node][word]

	for edge in bfs_Graph.edges():
		for word in ['weight']:
			bfs_Graph.edges[edge][word] = Graph.edges[edge][word]


	## -- plotting

	plt_graph_plots([bfs_Graph], partition_dict, color_dict, folder_path,'/bfs_graph.png')

	return bfs_Graph


def plt_graph_plots(Graph_list, partition_dict, color_dict, folder_path, file_name, alpha=.6, width_list = [3], color_list = ['black'], highlighted_nodes = None):

	fig, ax = plt.subplots(1, 1, figsize=(15, 15))

	#plotting the image (square by square)

	patches = []

	for key in partition_dict:

		square_edges = np.asarray([partition_dict[key][0]] + [partition_dict[key][2]] + [partition_dict[key][3]] + [
			partition_dict[key][1]] + [partition_dict[key][0]])

		s1 = Polygon(square_edges)
		patches.append(s1)

	p = PatchCollection(patches, alpha=alpha, cmap='YlOrRd', linewidth=.1, edgecolor='b')

	# getting colors
	colors = np.array(list(color_dict.values()))
	p.set_array(colors)
	ax.add_collection(p)

	#plotting graphs
	if len(Graph_list)!=len(width_list):
		width_list = len(Graph_list) * width_list
	if len(Graph_list)!=len(color_list):
		width_list = len(Graph_list) * color_list

	for Graph, width, color in zip(Graph_list,width_list,color_list):

		pos = nx.get_node_attributes(Graph, 'pos')
		nx.draw_networkx(Graph, pos, node_size=5, width=width, with_labels=False, edge_color=color,
						 alpha=0.5, ax=ax)

	if highlighted_nodes is not None:
		if len(highlighted_nodes)==1:

			nodes = highlighted_nodes[0]
			for i in range(len(nodes)):
				node = nodes[i]
				color = 'red'
				size = 4

				x = Graph_list[0].nodes[node]['pos'][0]
				y = Graph_list[0].nodes[node]['pos'][1]
				circle1 = plt.Circle((x, y), .01, color=color, fill=False, lw=size)
				ax.add_artist(circle1)
				# label = ax.annotate(str(node), xy=(x, y), fontsize=12, ha="center")


		elif len(highlighted_nodes)==2:

			nodes = highlighted_nodes[0]
			for i in range(len(nodes)):
				node = nodes[i]
				color = 'red'
				size = 4

				x = Graph_list[0].nodes[node]['pos'][0]
				y = Graph_list[0].nodes[node]['pos'][1]
				circle1 = plt.Circle((x, y), .01, color=color, fill=False, lw=size)
				ax.add_artist(circle1)
			# label = ax.annotate(str(node), xy=(x, y), fontsize=12, ha="center")

			nodes = highlighted_nodes[1]
			for i in range(len(nodes)):
				node = nodes[i]
				color = 'black'
				size = 2

				x = Graph_list[0].nodes[node]['pos'][0]
				y = Graph_list[0].nodes[node]['pos'][1]
				circle1 = plt.Circle((x, y), .01, color=color, fill=False, lw=size)
				ax.add_artist(circle1)
				# label = ax.annotate(str(node), xy=(x, y), fontsize=12, ha="center")

		else:
			print('more than two node lists to highlight!')

	plt.savefig(folder_path + file_name)
	plt.close('all')







def bfs_preprocess(image_path, new_size,number_of_colors, t1,t2, number_of_cc,graph_type,reversed_colors=True):
	'''
	This is the combination of pre_extraction_from_image and tree_approximation.
	:param image_path: string.
	:param number_of_colors: number of colors for the output image.
	:param t1: noise threshold. If t=0, then all the new pixels are preserved with their new colors.
	:param t2: threshold for the weights of the edges after pre-extraction.
	:param number_of_cc: number of connected components of the graph represented by the image. If None, then only 1
	cc is assumed.
	:param graph_type: 1 (to use edges and vertices of the grid), 2 (to use only edges).
	:return:
		bfs_Graph: bfs approximation of G.
	'''

	width, color_dict, folder_path = resizing_image(image_path, number_of_colors, new_size, t1, reversed_colors)

	Graph,_,_,_ = pre_extraction_from_image(image_path,new_size,t2,number_of_colors,number_of_cc,graph_type, t1, reversed_colors)

	partition_dict, _, _ = quality_measure.partition_set(width + 1)

	bfs_Graph = tree_approximation(Graph, color_dict, partition_dict, folder_path)

	with open(folder_path+'/bfs_extracted_graph.pkl', 'wb') as file:
		pkl.dump(bfs_Graph, file)

	with open(folder_path+'/real_part_dict.pkl', 'wb') as file:
		pkl.dump(partition_dict, file)

	return bfs_Graph



"""
Test
"""
#--------------------------- test 1 --------------------------------

new_size=100
ratio=new_size/1200
#print('ratio:',ratio)
t1=0.05
t2=.12
image_path = "./runs/graph_from_image/crop.png"
number_of_cc=1
number_of_colors=100
graph_type='1'

#Graph,_ = pre_extraction_from_image(image_path,new_size,t2,number_of_colors,number_of_cc,graph_type,t1)

#bfs_preprocess(image_path,new_size, number_of_colors, t1,t2, number_of_cc,graph_type )

#--------------------------- test 2 --------------------------------

new_size=100
ratio=new_size/1200
#print('ratio:',ratio)
t1=0.05
t2=.12
image_path = "./runs/graph_from_image/resized_crop.jpg"
number_of_cc=1
number_of_colors=100
graph_type='1'

#Graph,_ = pre_extraction_from_image(image_path,graph_type,t1,t2,number_of_colors,number_of_cc)

#bfs_preprocess(image_path, new_size,number_of_colors, t1,t2, number_of_cc,graph_type )

#----------------------------test 3---------------------------------------------------
new_size=100
ratio=new_size/1200
#print('ratio:',ratio)
t1=0.1
t2=.5
image_path = "./runs/graph_from_image/image_2.jpg"
number_of_cc=1
number_of_colors=100
graph_type='1'

#Graph,_ = pre_extraction_from_image(image_path,new_size,t2,number_of_colors,number_of_cc,graph_type,t1)

#bfs_preprocess(image_path,new_size, number_of_colors, t1,t2, number_of_cc,graph_type )


