import networkx as nx
from scipy.spatial import distance
import os
import random
import pickle as pkl
import numpy as np
import time
# -----------------------
import pre_extraction
import filtering
import quality_measure
#------------------------
nextrout_path = '/home/dtheuerkauf/Nextrout/python_scripts/'#open("./nextrout_path.txt", "r").readline()
print(nextrout_path)
#nextrout_path = '/is/sg2/dtheuerkauf/Nextrout/python_scripts/' #cluster path
def superimposing_graphs(graph_list, G_aux, union_type = None):

	unionG = nx.Graph()
	union_nodes = []
	union_edges = []

	for graph in graph_list:
		nodes = list(graph.nodes())
		edges = list(graph.edges())
		union_nodes = union_nodes + nodes
		union_edges = union_edges + edges

	union_nodes = list(set(union_nodes))
	union_edges = list(set(union_edges))
	unionG.add_nodes_from(union_nodes)
	unionG.add_edges_from(union_edges)

	for node in union_nodes:
		unionG.nodes[node]['pos'] = G_aux.nodes[node]['pos']

	for edge in union_edges:
		w_tdens = 0
		w_flux = 0
		avg_counter=0
		weight_flag = False # this is used to access either one or two weights when adding the edges
		for graph in graph_list:
			try:# if no error, then the edge *is* in the graph
				membership_test =graph.edges[edge]
				try: 
					w_tdens += graph.edges[edge[0], edge[1]]['tdens']
					
					avg_counter+=1
				except:
					w_tdens+= graph.edges[edge[0], edge[1]]['weight']
					avg_counter+=1
					weight_flag = True
			except:
				pass
		if weight_flag:
			unionG.edges[edge[0], edge[1]]['weight'] = w_tdens/avg_counter
		else:	
			unionG.edges[edge[0], edge[1]]['tdens'] = w_tdens/avg_counter
			
		unionG.edges[edge[0], edge[1]]['length'] = distance.euclidean(
			G_aux.nodes[edge[0]]['pos'], G_aux.nodes[edge[1]]['pos'])

		#unionG.edges[edge[0], edge[1]]['inv_weight'] = len(graph_list) / w
		#unionG.edges[edge[0], edge[1]]['iw_l'] = unionG.edges[edge[0], edge[1]]['length'] * unionG.edges[edge[0], edge[1]]['inv_weight']

	return unionG

def time_weights(G, beta_d):

	delta = (2-beta_d)/(3-beta_d)
	for edge in G.edges():
		l = G.edges[edge]['length']
		q = abs(G.edges[edge]['flux'])
		G.edges[edge]['time'] = l/(q**delta)

	return G

def image2net(image_path, N_runs, t2, t3, new_size = 'automatic',union_type = None, reversed_colors = True, noise_on=True, ds_interpolation = 'nearest', weighting_method_simplification = 'BPW', rseed = None):

	max_tol=5
	
	if rseed is None:

		rseed = list(range(max_tol))


	if len(rseed) < N_runs:

		print('not enough seeds')

	else:	

		
		# STEP 1: pre extraction

		# parameters:
		#t2 =.5 #<-- min
		#t3 =.6 #<--max
		number_of_colors = 0
		number_of_cc = 1
		graph_type = '1'
		t1 = 0.0 #threshold to clean the image (no real impact)

		#weighting_method_simplification = 'ER'

		print('step 1: computing Gpe.')

		G_pre_extracted, color_dict, partition_dict, folder_path = pre_extraction.pre_extraction_from_image(image_path,
																											new_size, t2,
																											number_of_colors,
																											number_of_cc,
																											graph_type, t1,
																											reversed_colors, t3,
																											ds_interpolation)
		
		print('connected',nx.is_connected(G_pre_extracted))

		# G_pre_extracted = filtering.bifurcation_paths(G,[])
		# adding some noise

		if noise_on:
			for edge in G_pre_extracted.edges():
				w = G_pre_extracted.edges[edge]['weight']
				G_pre_extracted.edges[edge]['weight']= w *random.uniform(0.99,1.01)

		
		print('graph stored at:'+folder_path)

		## STEP 2: tree approximation

		#parameters:
		beta_d = 1.5

		terminal_list = list(G_pre_extracted.nodes())

		#current_path = os.getcwd()
		#os.chdir(nextrout_path)

		print('step 2: computing Gtree.')

		G_tree = filtering.filtering_from_image(G_pre_extracted,
			                                           sources= terminal_list[0:1],
			                                           sinks = terminal_list[1:],
													   beta_d = beta_d,
													   color_dict = color_dict,
													   partition_dict = partition_dict,
													   weighting_method_simplification=weighting_method_simplification,
													   folder_name = folder_path)

		
		print('     is Gtree a tree?', str(nx.is_tree(G_tree)))

		G_tree = quality_measure.relabeling(G_tree, G_pre_extracted) # preprocessing needed for the filtering input format

		with open(folder_path + '/G_tree.pkl', 'wb') as file:
			pkl.dump(G_tree, file)

		#STEP #3: filtering

		

		# computing the leaves
		deg = nx.degree_centrality(G_tree)

		N = len(G_tree.nodes)

		for node in deg.keys():
			deg[node] = round(deg[node] * (N - 1))

		terminal_list = [node for node in G_tree.nodes() if deg[node] == 1]

		print('step 3: computing Gfs.')
		print('number of leaves:',len(terminal_list))


		Gf = {}
		sources = {}
		i=0
		max_=0
		while i<N_runs and max_<max_tol:

			rng = np.random.RandomState(seed=rseed[max_])

			print('i=',i,'/',N_runs)
			print('max_=',max_,'/',30)
			

			#getting a random index for the source

			random_source = rng.choice(terminal_list)
			print('chosen source:',random_source)
			rs_index = terminal_list.index(random_source)

			sources_rs = [terminal_list[entry] for entry in [rs_index]]
			sinks_rs = [node for node in terminal_list if node not in sources_rs]
			G_filtered = filtering.filtering_from_image(G_pre_extracted,
			                                           sources= sources_rs,
			                                           sinks = sinks_rs,
													   beta_d = beta_d,
													   color_dict = color_dict,
													   partition_dict = partition_dict,
													   weighting_method_simplification=weighting_method_simplification,
													   folder_name = folder_path)


			
			G_filtered = quality_measure.relabeling(G_filtered, G_pre_extracted)
			print('     is Gtree a tree?', str(nx.is_tree(G_filtered)))
			Gf[i] = G_filtered
			sources[i]=random_source
			i+=1
			max_+=1


		#os.chdir(current_path)

		pre_extraction.plt_graph_plots([G_tree],
										   partition_dict,
										   color_dict,
										   folder_path,
										   '/G_tree.png',
										   alpha=.6,
										   width_list=[3],
										   color_list=['black'])

		# superimposing the filtrations

		graph_list = list(Gf.values())

		for i in range(len(graph_list)):

			graph = graph_list[i]
			random_source = sources[i]
			print('------random source:'+str(random_source)+'----------: n_nodes:',len(graph.nodes()),', n_edges:',len(graph.edges()))


			pre_extraction.plt_graph_plots([graph],
										   partition_dict,
										   color_dict,
										   folder_path,
										   '/G_filtered' + str(random_source) + '.png',
										   alpha=.6,
										   width_list=[3],
										   color_list=['black'])
			with open(folder_path + '/G_filtered' + str(random_source) + '.pkl', 'wb') as file:
				pkl.dump(graph, file)



		G_filtered = superimposing_graphs(graph_list, G_pre_extracted, union_type)

		#adding value property

		#G_filtered = time_weights(G_filtered, beta_d)

		weights = [G_filtered.edges[edge]['tdens'] for edge in G_filtered.edges()]

		pre_extraction.plt_graph_plots([G_filtered],
									   partition_dict,
									   color_dict,
									   folder_path,
									   '/G_filtered.png',
									   alpha=.5,
									   width_list=[weights],
									   color_list=['black'],
									   highlighted_nodes = [ list(sources.values()) ])
		with open(folder_path + '/G_filtered.pkl', 'wb') as file:
			pkl.dump(G_filtered, file)

		#print(list(G_filtered.edges(data=True))[:10])

		parameters = [image_path, N_runs, new_size, t2, number_of_colors, number_of_cc, graph_type, t1, beta_d, t3, N_runs, max_, weighting_method_simplification, ds_interpolation]

		with open(folder_path + '/parameters.pkl', 'wb') as file:
			pkl.dump(parameters, file)
		#print(folder_path)
		#os.system('cp -r '+folder_path+' ../../data/output/test/'+folder_path.split('/')[-1]+'/nextrout/')
		try:
		  os.system('rm -r '+ '../../data/output/test/'+folder_path.split('/')[-1]+'/nextrout/')
		except:
		  pass

		#print(folder_path)

		#move_files(folder_path,'../../data/output/test/'+folder_path.split('/')[-1]+'/nextrout/')

		#os.system('cp -r '+folder_path + '  ' +folder_path+'nextrout_'+weighting_method_simplification)

		#os.system('rm -r '+folder_path)

		#print(G_filtered.edges(data=True))
		
		return G_filtered

def move_files(folder_origin, folder_dest):
	try:
		os.system('mkdir '+folder_dest)
	except:
		pass
	files = [folder_origin+'/'+f for f in os.listdir(folder_origin)]
	for f in files: 
		os.system('cp '+f+' '+folder_dest)


