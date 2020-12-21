import networkx as nx
import pickle as pkl
import random
from scipy.spatial import distance
from networkx.algorithms.approximation import steiner_tree
from networkx.algorithms import tree
import os
import time
import numpy as np
# ----------------------
import sys
import pre_extraction
import img2net
import filtering
import quality_measure
# -----------------------



def mst(image_path, N_runs, new_size = 'automatic',union_type = None, reversed_colors = True, noise_on=True, rseed = None, t2 =.1, ds_interpolation = 'nearest'):

	#cluster inputs
	cluster_flag = True
	G_filtered = 0

	if rseed is None:

		rseed = list(range(N_runs))

	if len(rseed)!= N_runs:

		print('not enough seeds')

	else:

		s = time.process_time()

		# STEP 1: pre extraction

		# parameters:
		#t2 = .5
		number_of_colors = 0
		number_of_cc = 1
		graph_type = '1'
		t1 = 0.0

		print(ds_interpolation)

		G_pre_extracted, color_dict, partition_dict, folder_path = pre_extraction.pre_extraction_from_image(image_path,
																											new_size, t2,
																											number_of_colors,
																											number_of_cc,
																											graph_type, t1,
																											reversed_colors,
																											ds_interpolation = ds_interpolation)
		Gpe_nodes = len(G_pre_extracted.nodes())
		Gpe_edges = len(G_pre_extracted.edges())
		print('Gpe: n_nodes,n_edges:',len(G_pre_extracted.nodes()),len(G_pre_extracted.edges()))

		#G_pre_extracted = filtering.bifurcation_paths(G_pre_extracted,[])

		## STEP 2: tree approximation
		
		#sys.exit()
		
		# bfs_Graph = steiner_tree(G_pre_extracted, list(G_pre_extracted.nodes()))
		mst = tree.maximum_spanning_edges(G_pre_extracted, algorithm='kruskal', data='weight')
		bfs_Graph = nx.Graph()
		bfs_Graph.add_edges_from(list(mst))

		# bfs_Graph = pre_extraction.tree_approximation(G_pre_extracted, color_dict, partition_dict, folder_path)

		with open(folder_path + '/bfs_extracted_graph.pkl', 'wb') as file:
			pkl.dump(bfs_Graph, file)

		with open(folder_path + '/other_par.pkl', 'wb') as file:
			pkl.dump([color_dict, partition_dict, folder_path], file)

		print('info stored')


		if cluster_flag:
			
			folder_path_cluster = './runs/'+image_path.split('/')[-2]+'/'+(image_path.split('/')[-1]).split('.')[0]
			'''
			print('fpc',folder_path_cluster)
			with open(folder_path_cluster + '/other_par.pkl', 'rb') as file:
				[color_dict, partition_dict, folder_path] = pkl.load(file)

			with open(folder_path_cluster + '/extracted_graph.pkl', 'rb') as file:
				G_pre_extracted = pkl.load(file)

			with open(folder_path_cluster + '/bfs_extracted_graph.pkl', 'rb') as file:
				bfs_Graph = pkl.load(file)
			'''
			# STEP #3: filtering (via MST)

			# computing the leaves
			deg = nx.degree_centrality(bfs_Graph)

			N = len(bfs_Graph.nodes)

			for node in deg.keys():
				deg[node] = round(deg[node] * (N - 1))

			terminal_candidates= [node for node in bfs_Graph.nodes() if deg[node] == 1]

			print('leaves',len(terminal_candidates))

			Gf = {}
			terminals = {}

			

			for i in range(N_runs):

				print('steiner tree n=',i)
				
				# getting a random terminals
				
				#terminal_list = [random.choice(terminal_candidates) for i in range(int(.025 * len(terminal_candidates)))]
				
				terminal_list = []

				print('random seed',rseed[i])

				rng = np.random.RandomState(seed=rseed[i])

				node = rng.choice(terminal_candidates)
				
				terminal_list.append(node)

				print('number of terminals',int(.025 * len(terminal_candidates)))
				
				while len(terminal_list)<int(.025 * len(terminal_candidates)):

					#print('i',i)
					
					node = rng.choice(terminal_candidates)
					#terminal_list.append(node)
					#print('node',node)
					node_pos = G_pre_extracted.nodes[node]['pos']
					dist_bool = []
					for elem in terminal_list:
						if elem != node:
							#print('elem!= node')
							elem_pos = G_pre_extracted.nodes[elem]['pos']
							val = distance.euclidean(elem_pos,node_pos)>.2
							dist_bool.append(val)
							#print('db',dist_bool)
					#print('db',dist_bool)
					if sum(dist_bool)==len(dist_bool):
						terminal_list.append(node)
				
				print('terminal_list',terminal_list)
				
				print('computing steiner tree:')
				G_filtered = steiner_tree(G_pre_extracted, terminal_list)
				
				G_filtered = quality_measure.relabeling(G_filtered, G_pre_extracted)

				Gf[i] = G_filtered
				
				terminals[i] = terminal_list
				
				#print(terminal_list)


			
			for i in range(N_runs):

				graph = Gf[i]
				terminal_list = terminals[i]
				#print('terminals',terminal_list)

				print('------random n=:'+str(i)+'----------: n_nodes:',len(graph.nodes()),', n_edges:',len(graph.edges()))


				pre_extraction.plt_graph_plots([graph],
											   partition_dict,
											   color_dict,
											   folder_path_cluster,
											   '/G_filtered' + str(i) + '.png',
											   alpha=.6,
											   width_list=[3],
											   color_list=['black'],
											   highlighted_nodes=[terminal_list])
				#with open(folder_path + '/G_filtered' + str(i) + '.pkl', 'wb') as file:
				#    pkl.dump(graph, file)

			G_filtered = img2net.superimposing_graphs(Gf.values(), G_pre_extracted, union_type)

			weights = [8*G_filtered.edges[edge]['weight'] for edge in G_filtered.edges()]

			pre_extraction.plt_graph_plots([G_filtered],
										   partition_dict,
										   color_dict,
										   folder_path,
										   '/G_filtered.png',
										   alpha=.5,
										   width_list=[weights],
										   color_list=['black']
										   )
			with open(folder_path_cluster + '/G_filtered.pkl', 'wb') as file:
				pkl.dump(G_filtered, file)
			
			print('storing results at:',folder_path_cluster)
			'''
			try:
				os.system('rm -r '+ '../../data/output/test/'+folder_path.split('/')[-1]+'/mst/')
			except:
				pass
			'''
			#img2net.move_files(folder_path_cluster,'../../data/output/test/'+folder_path.split('/')[-1]+'/mst/')

			#os.system('rm -r '+folder_path)

			e = time.process_time()

			execution_time = e-s
			'''
			with open( '../../data/output/test/'+folder_path.split('/')[-1]+'/mst/ex_time.pkl', 'wb') as file:
				pkl.dump(execution_time, file)
			'''
		return G_filtered, Gpe_nodes, Gpe_edges



for test_flag in ['others']:#,1,2,3,5]:

	if test_flag == 1:

		times = []

		t2 =.45

		path = './test_images/crop.png'


		mst(path,N_runs = 5,new_size=180,  t2 = t2)
 

	elif test_flag == 2:

		t2 = .24

		path = './test_images/084_ro0117_f7-2.jpg'

		mst(path, N_runs=5, new_size=100, reversed_colors= True, t2 = t2)

	elif test_flag == 3:

		t2 = .19
		path = './test_images/venezuelan_rivers.png'

		mst(path, N_runs=5, new_size=170, reversed_colors= True,  t2 =t2)

	elif test_flag == 4:

		t2 = .1
		path = './test_images/eye_vein.jpg'

		mst(path, N_runs=5, new_size=200, reversed_colors= True,  t2 = t2)

	elif test_flag == 5:

		t2 = .45
		path = './test_images/image_1626.jpg'

		mst(path, N_runs=5, new_size=150, reversed_colors= True,  t2 = t2)

	elif test_flag == 'retinal':

		imgs = [
		'../../data/input/retinal_vessel/stare-images/'+folder
		for folder in os.listdir('../../data/input/retinal_vessel/stare-images/') if 'DS_Store' not in folder and  'crop' in folder
	]


		for image_path in imgs:
			print(os.getcwd())
			mst(image_path, N_runs=5, new_size = 200, reversed_colors =True, ds_interpolation = 'nearest')

	elif test_flag == 'others':
		
		
		'''
		times = []
		imgs = [
		'../../data/input/physarum/'+folder
		for folder in os.listdir('../../data/input/physarum/') if 'DS_Store' not in folder and 'crop' in folder 
		]	 
		
		for image_path in imgs[:1]:
			print(image_path)
			start = time.time()
			mst(image_path, N_runs=5, new_size = 300, reversed_colors =True, ds_interpolation = 'nearest', t2 = .5)
			end = time.time()-start
			times.append(end)
		with open('mst_times.pkl', 'wb') as file:
			pkl.dump(times, file)
		'''
		
		imgs = [
		'../../data/input/rivers/'+folder
		for folder in os.listdir('../../data/input/rivers/') 
		if 'DS_Store' not in folder 
		and 'crop' in folder
		and ('zambian' in folder
		or 'papua' in folder
		or 'angolan' in folder
		or 'canadian' in folder
		or 'yenisei' in folder)
		]	
		#times = []
		
		print(imgs)
		
		for image_path in imgs:
			
			print(image_path)
			
			start = time.time()
			_,nodes, edges = mst(image_path, N_runs=5, new_size = 200, reversed_colors =True, ds_interpolation = 'area', t2 = .2)
			ArithmeticError(args)end = time.time()
			t = end-start
			
			with open('mst_times'+image_path.split('/')[-1].split('.')[0]+'.pkl', 'wb') as file:
				pkl.dump([image_path,t,nodes, edges], file)
		

