import os
import networkx as nx
import numpy as np
import math
import itertools
from scipy.spatial import distance
import matplotlib.pyplot as plt
###-------------
import utils
###-------------



def source_rect_cnst_test( x, y):
    x1, y1, x2, y2 = 1.0 / 8.0, 1 / 4, 3 / 8, 3 / 4
    if (x > x1 and x < x2 and y > y1 and y < y2):
        return True
    else :
        return False


def sink_rect_cnst_test( x, y):
    x3, y3, x4, y4 = 5.0 / 8.0, 1 / 4, 7 / 8, 3 / 4
    if (x > x3 and x < x4 and y > y3 and y < y4):
        return True
    else :
        return False


def fsource(x,y,flag):
    '''
    Source selection based on coordinates of a node. If f is the source/sink function and z=f(x,y)!=0, then (x,y)
    is source.
    :param x: x coordinate.
    :param y: y coordinate.
    :param flag: string for source: '3rc',
                                    '2rcs',
                                    '2rcl',
                                    '3rch',
                                    '3rcl',
                                    '5rch',
                                    '5rcm',
                                    '6rcl',
                                    '5rcl',
                                    '4rch',
                                    '4rcl',
                                    '4rcm',...
    :return:
        z: source value.
    '''
    radio=.05
    z=0
    #------------------------------------------------------------------
    if flag == '1_4_small_squares':
        if x<.1 and y<.1:
            z=1
    #------------------------------------------------------------------
    if flag =='signal_prop':
        # taken from https://www.pnas.org/content/114/20/5136
        if (x > .2 and x < .8 and y > .5 and y < .7):
            z=1
    #-----------------------------------------------------------------
    if flag == '1_rect' or flag == '3_rect' or flag == '1_obstacles' or flag == '3_obstacles':
        if math.sqrt((x-.2)**2+(y-.5)**2)<=radio:
            z=1
        if flag == '1_obstacles' or flag == '3_obstacles':
            if math.sqrt((x-.5)**2+(y-.5)**2)<=radio:
                z=1
    #------------------------------------------------------------------------------------
    if flag=='2_rect' or flag == '3_rect' or flag == '2_obstacles' or flag == '3_obstacles':
        if math.sqrt((x-.2)**2+(y-.2)**2)<=radio or math.sqrt((x-.2)**2+(y-.8)**2)<=radio:
            z=1
        if flag == '2_obstacles' or flag == '3_obstacles':
            if math.sqrt((x-.4)**2+(y-.3)**2)<=radio or math.sqrt((x-.4)**2+(y-.7)**2)<=radio:
                z=1    
    #-----------------------------------------------------------------------------------
    #a circle in the center of the domain
    if flag =='center':

        radio=.01

        if math.sqrt((x-.5)**2+(y-.5)**2)<=radio:
            z=10 
    if flag == 'rc':
        if 0.9<=y:
            z=1

    #-----------------------------------------------------------------------------------            
    elif flag =='top_circle':
        radio = .01
        if math.sqrt((x-.5)**2+(y-.8)**2)<=.05:
            z=1
    #-----------------------------------------------------------------------------------            
    elif flag =='pp_paper':
        radio = .05
        if (math.sqrt((x-.5)**2+(y-.8)**2)<=2*radio 
        or math.sqrt((x-.2)**2+(y-.8)**2)<=2*radio
        or math.sqrt((x-.8)**2+(y-.8)**2)<=2*radio 
        or math.sqrt((x-.0)**2+(y-.6)**2)<=2*radio 
        or math.sqrt((x-.0)**2+(y-1)**2)<=2*radio
        or math.sqrt((x-.8)**2+(y-1)**2)<=2*radio  
        or math.sqrt((x-.4)**2+(y-1)**2)<=2*radio 
        or math.sqrt((x-.7)**2+(y-.9)**2)<=2*radio 
        or math.sqrt((x-.35)**2+(y-.55)**2)<=2*radio 
        or math.sqrt((x-.35)**2+(y-1)**2)<=2*radio
        or math.sqrt((x-1)**2+(y-1)**2)<=3*radio):
            z=1


    #-----------------------------------------------------------------------------------
    # 3 rectangles, 2 as sources and 1 as sink placed in the middle
    x1 , y1 , x2 , y2, y3, y4 = 0.125, 0.125, 0.375, 0.4, 0.650, 0.9
    if flag=='3rc':
        if ((x > x1 and x < x2 and y > y1 and y < y2) or (x > x1 and x < x2 and y > y3 and y < y4)):
            z=1
    #-----------------------------------------------------------------------------------
            
    ### 2 rectangles, 1 as source, 1 as sink placed at bottom, top
    x5 , y5 , x6 , y6 = 0.125, 0.125, 0.375, 0.4
    if flag=='2rcs':
        if ((x > x5 and x < x6 and y > y5 and y < y6)):
            z=1
    #-----------------------------------------------------------------------------------
            
    ### 2 rectangles, 1 as source, 1 as sink placed at bottom, top
    x7 , y7 , x8 , y8 = 0.125, 0.125, 0.375, 0.4
    if flag=='2rcl':
        if ((x > x7 and x < x8 and y > y7 and y < y8)):
            z=1
    #-----------------------------------------------------------------------------------
            
    ### 3 rectangles, 2 as sources, 1 sink placed in the higher part
    x9 , y9 , x10 , y10, y11, y12 = 0.125, 0.125, 0.375, 0.4, 0.650, 0.9
    if flag=='3rch':
        if ((x > x9 and x < x10 and y > y9 and y < y10) or (x > x9 and x < x10 and y > y11 and y < y12)):
            z=1
    #-----------------------------------------------------------------------------------
            
    ### 3 rectangles, 2 as sources, 1 sink placed in the lower part
    x13 , y13 , x14 , y14, y15, y16 = 0.125, 0.125, 0.375, 0.4, 0.650, 0.9
    if flag=='3rcl':
        if ((x > x13 and x < x14 and y > y13 and y < y14) or (x > x13 and x < x14 and y > y15 and y < y16)):
            z=1
   

    #-----------------------------------------------------------------------------------
    x17 ,  x18 , y17, y18, y19, y20, y21, y22, x19, x20= 0.125, 0.225, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825, 0.45, 0.65
    if flag=='5rcm':
        if ((x > x17 and x < x18 and y > y17 and y < y18) or (x > x17 and x < x18 and y > y19 and y < y20) or (x > x17 and x < x18 and y > y21 and y < y22) or (x > x19 and x < x20 and y > y19 and y < y20)):
            z=1

    #-----------------------------------------------------------------------------------
    x21 ,  x22 , y17, y18, y19, y20, y21, y22, x19, x20= 0.125, 0.225, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825, 0.45, 0.65
    if flag=='5rch':
        if ((x > x21 and x < x22 and y > y17 and y < y18) or (x > x17 and x < x18 and y > y19 and y < y20) or (x > x17 and x < x18 and y > y21 and y < y22) or (x > x19 and x < x20 and y > y19 and y < y20)):
            z=1
    #-----------------------------------------------------------------------------------
    x21 ,  x22 , y17, y18, y19, y20, y21, y22, x19, x20= 0.05, 0.15, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825, 0.825, 0.95
    if flag=='5rch2':
        if ((x > x21 and x < x22 and y > y17 and y < y18) or (x > x21 and x < x22 and y > y19 and y < y20) or (x > x21 and x < x22 and y > y21 and y < y22) or (x > x19 and x < x20 and y > y19 and y < y20)):
            z=1
    #-----------------------------------------------------------------------------------
    
    x21 ,  x22 , y17, y18, y19, y20, y21, y22, x19, x20= 0.125, 0.225, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825, 0.45, 0.65
    if flag=='5rcl':
        if ((x > x21 and x < x22 and y > y17 and y < y18) or (x > x17 and x < x18 and y > y19 and y < y20) or (x > x17 and x < x18 and y > y21 and y < y22) or (x > x19 and x < x20 and y > y19 and y < y20)):
            z=1

 
    #-----------------------------------------------------------------------------------
    x17 ,  x18 , y17, y18, y19, y20, y21, y22= 0.125, 0.225, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825
    if flag=='4rcl':
        if ((x > x17 and x < x18 and y > y17 and y < y18) or (x > x17 and x < x18 and y > y19 and y < y20) or (x > x17 and x < x18 and y > y21 and y < y22)):
            z=1

     #-----------------------------------------------------------------------------------
    x17 ,  x18 , y17, y18, y19, y20, y21, y22, = 0.125, 0.225, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825
    if flag=='4rcm':
        if ((x > x17 and x < x18 and y > y17 and y < y18) or (x > x17 and x < x18 and y > y19 and y < y20) or (x > x17 and x < x18 and y > y21 and y < y22)):
            z=1

     #-----------------------------------------------------------------------------------
    x17 ,  x18 , y17, y18, y19, y20, y21, y22= 0.125, 0.225, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825
    if flag=='4rch':
        if ((x > x17 and x < x18 and y > y17 and y < y18) or (x > x17 and x < x18 and y > y19 and y < y20) or (x > x17 and x < x18 and y > y21 and y < y22)):
            z=1

 
    #-----------------------------------------------------------------------------------
    return z



def fsink(x,y,flag):
    '''
    Sink selection based on coordinates of a node. If f is the source/sink function and z=f(x,y)!=0, then (x,y)
    is sink.
    :param x: x coordinate.
    :param y: y coordinate.
    :param flag: string for source: '3rc',
                                    '2rcs',
                                    '2rcl',
                                    '3rch',
                                    '3rcl',
                                    '5rch',
                                    '5rcm',
                                    '6rcl',
                                    '5rcl',
                                    '4rch',
                                    '4rcl',
                                    '4rcm', ...
    :return:
        z: sink value.
    '''
    radio=.05
    z=0
    #----------------------------------------------------------------------------------
    if flag == '1_4_small_squares':
        if (( x < .1 and y > .9 ) or  #left-top
         (x > .9 and y > .9 ) or  #right-top
        ( x > .9 and y < .1) or  #right-bottom
        (x > .45 and x < .55 and y > .45 and y < .55) ): #center
            z=1
    #----------------------------------------------------------------------------------
    if flag =='4_circles':
        if (math.sqrt((x-.35)**2+(y-.35)**2)<=radio) or (math.sqrt((x-.35)**2+(y-.65)**2)<=radio) or (math.sqrt((x-.65)**2+(y-.35)**2)<=radio) or (math.sqrt((x-.65)**2+(y-.65)**2)<=radio):
            z=1
    #-----------------------------------------------------------------------------------
    if flag =='signal_prop':
        # taken from https://www.pnas.org/content/114/20/5136
        if (x > .2 and x < .8 and y > .1 and y < .5):
            z=1

    #-----------------------------------------------------------------------------------
    if flag == '1_rect':
        if math.sqrt((x-.8)**2+(y-.5)**2)<=radio:
            z=2


    #-----------------------------------------------------------------------------------
    elif flag == '12cl':
        y0=0.04
        if (math.sqrt((x-.05)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.1)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.15)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.20)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.25)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.30)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.35)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.40)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.45)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.50)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.55)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.60)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.65)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.70)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.75)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.80)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.85)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.90)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.95)**2+(y-y0)**2)<=radio) or (math.sqrt((x-.99)**2+(y-y0)**2)<=radio):
            z=1            
    #-----------------------------------------------------------------------------------
    elif flag == 'rc':
        if not 0.9<=y:
            z=1
    #-----------------------------------------------------------------------------------
    elif flag == '2_lines':
        if 0.1>=y or .5<=y<=.6:
            z=1
    #-----------------------------------------------------------------------------------
    elif flag == '6_lines':
        if 0.05>=y or .25<=y<=.3 or .35<=y<=.4 or .45<=y<=.5 or .55<=y<=.6 or .65<=y<=.7:
            z=1  
    #-----------------------------------------------------------------------------------            
    elif flag =='low_circle':
        radio = .01
        if math.sqrt((x-.5)**2+(y-.2)**2)<=.05:
            z=1
    #-----------------------------------------------------------------------------------            
    elif flag =='pp_paper':
        radio = .05
        if (math.sqrt((x-.5)**2+(y-.2)**2)<=2*radio 
        or math.sqrt((x-.8)**2+(y-.2)**2)<=2*radio 
        or math.sqrt((x-.2)**2+(y-.2)**2)<=2*radio
        or math.sqrt((x-.35)**2+(y-.1)**2)<=2*radio
        or math.sqrt((x-.1)**2+(y-.1)**2)<=2*radio
        or math.sqrt((x-.75)**2+(y-.0)**2)<=2*radio
        or math.sqrt((x-1)**2+(y-.0)**2)<=6*radio):
            z=1

    #-----------------------------------------------------------------------------------            
    elif flag =='circle':
        radio = .01
        if math.sqrt((x-.5)**2+(y-.5)**2)<=.45 and not math.sqrt((x-.5)**2+(y-.5)**2) <=radio :
            z=1
    #--------------------------------------------------------------------------            

    x1 , x2 , y1, y2 = 0.625, 0.875, 0.45, 0.575

    if flag=='3rc':
        if (x > x1 and x < x2 and y > y1 and y < y2): 
            z=1
    #-----------------------------------------------------------------------------------            

    x3 , x4 , y3, y4 = 0.625, 0.875, 0.25, 0.38

    if flag=='3rcl':
        if (x > x3 and x < x4 and y > y3 and y < y4): 
            z=1
   #-----------------------------------------------------------------------------------            

    x5 , x6 , y5, y6 = 0.625, 0.875, 0.75, 0.9

    if flag=='3rch':
        if (x > x5 and x < x6 and y > y5 and y < y6): 
            z=1
    #-----------------------------------------------------------------------------------

    x7 , x8 , y7, y8 = 0.625, 0.875, 0.75, 0.9

    if flag=='2rcs':
        if (x > x7 and x < x8 and y > y7 and y < y8): 
            z=1

    #-----------------------------------------------------------------------------------        
    x9 , x10 , y9, y10 = 0.625, 0.875,  0.25, 0.45

    if flag=='2rcl':
        if (x > x9 and x < x10 and y > y9 and y < y10): 
            z=1
    #-----------------------------------------------------------------------------------        
    x11 , x12 , y11, y12 = 0.725, 0.925, 0.4, 0.7

    if flag=='5rcm':
        if (x > x11 and x < x12 and y > y11 and y < y12): 
            z=1

    #-----------------------------------------------------------------------------------        
    x13 , x14 , y13, y14 = 0.825, 0.95, 0.8, 0.95
    x01 , x02 , y01, y02 = 0.825, 0.95, 0.1, 0.225
    if flag=='6rcl':
        if (x > x13 and x < x14 and y > y13 and y < y14) or (x > x01 and x < x02 and y > y01 and y < y02): 
            z=1

    #-----------------------------------------------------------------------------------        
    x13 , x14 , y13, y14 = 0.725, 0.925, 0.75, 0.9

    if flag=='5rch':
        if (x > x13 and x < x14 and y > y13 and y < y14): 
            z=1
    #-----------------------------------------------------------------------------------        
    x15 , x16 , y15, y16 = 0.725, 0.925,  0.25, 0.38

    if flag=='5rcl':
        if (x > x15 and x < x16 and y > y15 and y < y16): 
            z=1
   
   #-----------------------------------------------------------------------------------     
    x17 , x18 , y17, y18 = 0.625, 0.875, 0.4, 0.7

    if flag=='4rcm':
        if (x > x17 and x < x18 and y > y17 and y < y18): 
            z=1

    #-----------------------------------------------------------------------------------        
    x13 , x14 , y13, y14 = 0.725, 0.925, 0.75, 0.9

    if flag=='4rch':
        if (x > x13 and x < x14 and y > y13 and y < y14): 
            z=1
    #-----------------------------------------------------------------------------------        
    x15 , x16 , y15, y16 = 0.625, 0.875,  0.25, 0.38

    if flag=='4rcl':
        if (x > x15 and x < x16 and y > y15 and y < y16): 
            z=1

    #-----------------------------------------------------------------------------------
    return z


def source_graph(X_bar,flag):
    nnzero=0
    source_values=[]
    for node in X_bar.nodes():
        x,y = X_bar.nodes[node]['pos']
        #print(x,y)
        z=fsource(x,y,flag)
        X_bar.nodes[node]['source']=z
        if z!=0:
            nnzero+=1
        source_values.append(z)
    return X_bar,source_values,nnzero

def sink_graph(X_bar,flag):
    nnzero=0
    sink_values=[]
    for node in X_bar.nodes():
        x,y = X_bar.nodes[node]['pos']
        #print(x,y)
        z=fsink(x,y,flag)
        X_bar.nodes[node]['sink']=z
        if z!=0:
            nnzero+=1
        sink_values.append(z)
    return X_bar,sink_values,nnzero

def source_sink_generator(folder_name, ndiv, flag_source, flag_sink):
    #based on jupyter notebook: Source_and_Sink_functions
    # requirements: graph_cell_folder
    #Reading the file
    print('Loading the graph cell information.')
    file1_ = open("./graph_cell_folder/graph_cell_"+str(ndiv)+".dat", "r")
    graph_coord_triang = file1_.readlines()
    file1_.close()
    n_nodes=int(graph_coord_triang[0][:12])
    graph_coordinates=graph_coord_triang[2:2+n_nodes]
    #Getting the positions of the nodes
    bar_pos = utils.bar2dict(graph_coordinates)
    #Defining the positions of the nodes
    X_bar=nx.Graph()
    X_bar.add_nodes_from(bar_pos.keys())
    for n, p in bar_pos.items():
        X_bar.nodes[n]['pos'] = p
    pos=nx.get_node_attributes(X_bar,'pos')
    
    if flag_source != 'rect_cnst':
        #Generating source graph
        print('Generating source graph.')
        X_source=X_bar.copy()
        _, source_values, nnzero_source = source_graph(X_source,flag_source)
        #Writing the source to a file
        print('Saving source into source.dat.')
        f= open(folder_name+"/input/source.dat","w+")
        n_nodes=len(X_bar.nodes.keys())
        print('number of nodes',n_nodes)
        f.write("1"+(7-len(str(n_nodes)))*" "+str(n_nodes)+"\n")
        f.write("time    0.0"+"\n")
        f.write(str(nnzero_source)+"\n")
        for node in X_source.nodes():
            z=X_source.nodes[node]['source']
            if z!=0:
                f.write("  "+node+" "+str(float(z))+"\n")
        f.write("time 1e30")
        f.close() 

    if flag_sink != 'rect_cnst':
        #Generating sink graph
        print('Generating sink graph.')
        X_sink=X_bar.copy()
        _, sink_values, nnzero_sink = sink_graph(X_sink,flag_sink)
        #Writing the sink to a file
        print('Saving sink into sink.dat')
        f= open(folder_name+"/input/sink.dat","w+")
        n_nodes=len(X_bar.nodes.keys())
        f.write("1"+(7-len(str(n_nodes)))*" "+str(n_nodes)+"\n")
        f.write("time    0.0"+"\n")
        f.write(str(nnzero_sink)+"\n")
        for node in X_sink.nodes():
            z=X_sink.nodes[node]['sink']
            if z!=0:
                f.write("  "+node+" "+str(float(z))+"\n")
        f.write("time 1e30")
        f.close() 


def source_sink_preprocess(folder_name):

    print('Source/sink preprocess.')
    sink = folder_name+'/input/sink.dat'
    source = folder_name+'/input/source.dat'
    new_forcing = folder_name+'/input/new_forcing.dat'
    command = '../otp_utilities/globals/axpy_timedata/axpy.out -1.0 ' + sink +'  '+source+'  '+new_forcing
    print('Computing the new forcing term.')
    os.system(command)
    subgrid = folder_name+'/input/subgrid.dat'
    parent = folder_name+'/input/parent.dat'
    new_forcing_subgrid = folder_name+'/input/new_forcing_subgrid.dat'
    command = '../otp_utilities/geometry/interpolate_timedata/interpolate_timedata.out'+'  '+subgrid+'  '+parent+'  '+new_forcing+'  '+new_forcing_subgrid
    print('Interpolating the new forcing onto the subgrid.')
    os.system(command)
    rhs = folder_name+'/input/rhs.dat'
    print('Generating the right hand side (rhs).')
    command = '../otp_utilities/p1galerkin/makerhs/makerhs.out'+' '+subgrid+' '+new_forcing_subgrid+' '+rhs
    os.system(command)
    rhs_integrated = folder_name+'/input/rhs_integrated.dat'
    print('Checking balance of the rhs.')
    command = 'python ../otp_utilities/globals/python_timedata/balance.py'+' '+rhs+'  '+rhs_integrated
    os.system(command)
    grid = folder_name+'/input/grid.dat'
    source_vtk = folder_name+'/input/source.vtk'
    sink_vtk = folder_name+'/input/sink.vtk'
    print('Generating source.vtk and sink.vtk.')
    command = '../otp_utilities/geometry/timedata2vtk/timedata2vtk.out'+' '+grid+'  '+source+' '+source_vtk
    os.system(command)
    command = '../otp_utilities/geometry/timedata2vtk/timedata2vtk.out'+' '+grid+'  '+sink+' '+sink_vtk
    os.system(command)
    


def source_plot(flag,ax=None):
    return_ax = True

    if flag == '1_4_small_squares':
        x1 = .0
        x2 = .1 
        y1 = .0
        y2 = .1
        _ , _, ax = rect_plot(x1,x2,y1,y2,ax,'source')
    elif flag =='signal_prop':
        x1 , y1 , x2 , y2 = .2,.5,.8, .7
        _ , _, ax = rect_plot(x1,x2,y1,y2,ax,'source')
    elif flag == 'rect_cnst':
        x1 , y1 , x2 , y2 = 1.0/8.0,1/4, 3/8, 3/4
        _ , _, ax = rect_plot(x1,x2,y1,y2,ax,'source')
    elif flag=='2rcs':
        x1 , y1 , x2 , y2 = 0.125, 0.125, 0.375, 0.4
        _ , _, ax = rect_plot(x1,x2,y1,y2,ax,'source')
    elif flag=='3rc':
        x1 , y1 , x2 , y2, y3, y4 = 0.125, 0.125, 0.375, 0.4, 0.650, 0.9
        _ , _, ax = rect_plot(x1,x2,y1,y2,ax,'source')
        _ , _, ax = rect_plot(x1,x2,y3,y4,ax,'source')
    elif flag=='2rcl':
        x7 , y7 , x8 , y8 = 0.125, 0.125, 0.375, 0.4
        _ , _, ax = rect_plot(x7,x8,y7,y8,ax,'source')
    elif flag=='3rch':
        x9 , y9 , x10 , y10, y11, y12 = 0.125, 0.125, 0.375, 0.4, 0.650, 0.9
        _ , _, ax = rect_plot(x9,x10,y9,y10,ax,'source')
        _ , _, ax = rect_plot(x9,x10,y11,y12,ax,'source')
    elif flag=='3rcl':
        x13 , y13 , x14 , y14, y15, y16 = 0.125, 0.125, 0.375, 0.4, 0.650, 0.9
        _ , _, ax = rect_plot(x13,x14,y13,y14,ax,'source')
        _ , _, ax = rect_plot(x13,x14,y15,y16,ax,'source')
    elif flag=='5rcm':
        x17 ,  x18 , y17, y18, y19, y20, y21, y22, x19, x20= 0.125, 0.225, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825, 0.45, 0.65
        _ , _, ax = rect_plot(x17,x18,y17,y18,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y19,y20,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y21,y22,ax,'source')
        _ , _, ax = rect_plot(x19,x20,y19,y20,ax,'source')
    elif flag=='5rch':  
        x17 ,  x18, x21 ,  x22 , y17, y18, y19, y20, y21, y22, x19, x20=0.125, 0.225,  0.125, 0.225, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825, 0.45, 0.65
        _ , _, ax = rect_plot(x21,x22,y17,y18,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y19,y20,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y21,y22,ax,'source')
        _ , _, ax = rect_plot(x19,x20,y19,y20,ax,'source')
    elif flag=='6rcl':
        x17 ,  x18, y17, y18, y19, y20, y21, y22, x19, x20= 0.05, 0.15, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825, 0.825, 0.95
        _ , _, ax = rect_plot(x17,x18,y17,y18,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y19,y20,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y21,y22,ax,'source')
        _ , _, ax = rect_plot(x19,x20,y19,y20,ax,'source')
    elif flag=='5rcl':
        x17 ,  x18, x21 ,  x22 , y17, y18, y19, y20, y21, y22, x19, x20= 0.125, 0.225, 0.125, 0.225, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825, 0.45, 0.65
        _ , _, ax = rect_plot(x21,x22,y17,y18,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y19,y20,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y21,y22,ax,'source')
        _ , _, ax = rect_plot(x19,x20,y19,y20,ax,'source')
    elif flag=='4rcl':
        x17 ,  x18 , y17, y18, y19, y20, y21, y22= 0.125, 0.225, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825
        _ , _, ax = rect_plot(x17,x18,y17,y18,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y19,y20,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y21,y22,ax,'source')
    elif flag=='4rcm':
        x17 ,  x18 , y17, y18, y19, y20, y21, y22, = 0.125, 0.225, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825
        _ , _, ax = rect_plot(x17,x18,y17,y18,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y19,y20,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y21,y22,ax,'source')
    elif flag=='4rch':
        x17 ,  x18 , y17, y18, y19, y20, y21, y22= 0.125, 0.225, 0.125, 0.3, 0.4, 0.575, 0.7, 0.825
        _ , _, ax = rect_plot(x17,x18,y17,y18,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y19,y20,ax,'source')
        _ , _, ax = rect_plot(x17,x18,y21,y22,ax,'source')
    elif flag == '1_rect' or flag == '3_rect'or flag == '1_obstacles' or flag == '3_obstacles':
        return_ax= True
        circle1=plt.Circle((.2, .5), 0.05, color='green', fill=False)
        ax.add_artist(circle1)
        if flag == '1_obstacles' or flag == '3_obstacles':
            circle2=plt.Circle((.5, .5), 0.05, color='green', fill=False)
            ax.add_artist(circle2)
    elif flag=='2_rect' or flag == '3_rect' or flag == '2_obstacles' or flag == '3_obstacles':
        return_ax= True
        circle11=plt.Circle((.2, .2), 0.05, color='green', fill=False)
        circle12=plt.Circle((.2, .8), 0.05, color='green', fill=False)
        ax.add_artist(circle11)
        ax.add_artist(circle12)
        if flag == '2_obstacles' or flag == '3_obstacles':
            circle21=plt.Circle((.4, .3), 0.05, color='green', fill=False)
            circle22=plt.Circle((.4, .7), 0.05, color='green', fill=False)
            ax.add_artist(circle21)
            ax.add_artist(circle22)
    elif flag =='center':
        return_ax= True  
        circle1=plt.Circle((.5, .5), .01, color='green', fill=False)
        ax.add_artist(circle1)
    else:
        print('Flag not defined!')
    if return_ax == True:
        return ax

def sink_plot(flag,ax=None):
    return_ax = True

    if flag == '1_4_small_squares':
        _ , _, ax = rect_plot(.0,.1,.9,1,ax,'sink') #x1,x2,y1,y2
        _ , _, ax = rect_plot(.9,1,.9,1,ax,'sink')
        _ , _, ax = rect_plot(.9,1,0,.1,ax,'sink')
        _ , _, ax = rect_plot(.45,.55,.45,.55,ax,'sink')
    elif flag =='signal_prop':
        x1 , y1 , x2 , y2 = .2,.1,.8, .5
        _ , _, ax = rect_plot(x1,x2,y1,y2,ax,'sink')
    elif flag == 'rect_cnst':
        x3 , y3 , x4 , y4 = 5.0/8.0,1/4, 7/8, 3/4
        _ , _, ax = rect_plot(x3,x4,y3,y4,ax,'sink')

    elif flag =='circle':
        return_ax= True    
        circle1=plt.Circle((.5, .5), 0.45, color='red', fill=False)
        ax.add_artist(circle1)
    elif flag=='3rc':
        x1 , x2 , y1, y2 = 0.625, 0.875, 0.45, 0.575

        _ , _, ax = rect_plot(x1,x2,y1,y2,ax,'sink')
    elif flag=='3rcl': 
        x3 , x4 , y3, y4 = 0.625, 0.875, 0.25, 0.38
        _ , _, ax = rect_plot(x3,x4,y3,y4,ax,'sink')
    elif flag=='3rch':
        x5 , x6 , y5, y6 = 0.625, 0.875, 0.75, 0.9
        _ , _, ax = rect_plot(x5,x6,y5,y6,ax,'sink')
    elif flag=='2rcs':
        x7 , x8 , y7, y8 = 0.625, 0.875, 0.75, 0.9
        _ , _, ax = rect_plot(x7,x8,y7,y8,ax,'sink')
    elif flag=='2rcl':
        x9 , x10 , y9, y10 = 0.625, 0.875,  0.25, 0.45
        _ , _, ax = rect_plot(x9,x10,y9,y10,ax,'sink')
    elif flag=='5rcm':
        x11 , x12 , y11, y12 = 0.725, 0.925, 0.4, 0.7
        _ , _, ax = rect_plot(x11,x12,y11,y12,ax,'sink')
    elif flag=='5rch':
        x13 , x14 , y13, y14 = 0.725, 0.925, 0.75, 0.9
        _ , _, ax = rect_plot(x13,x14,y13,y14,ax,'sink')
    elif flag=='6rcl':
        x13 , x14 , y13, y14 = 0.825, 0.95, 0.8, 0.95
        _ , _, ax = rect_plot(x13,x14,y13,y14,ax,'sink')
        x13 , x14 , y13, y14 = 0.825, 0.95, 0.1, 0.225,
        _ , _, ax = rect_plot(x13,x14,y13,y14,ax,'sink')
    elif flag=='5rcl':
        x15 , x16 , y15, y16 = 0.725, 0.925,  0.25, 0.38
        _ , _, ax = rect_plot(x15,x16,y15,y16,ax,'sink')
    elif flag=='4rcm':
        x17 , x18 , y17, y18 = 0.625, 0.875, 0.4, 0.7
        _ , _, ax = rect_plot(x17,x18,y17,y18,ax,'sink')
    elif flag=='4rch':
        x13 , x14 , y13, y14 = 0.725, 0.925, 0.75, 0.9
        _ , _, ax = rect_plot(x13,x14,y13,y14,ax,'sink')            
    elif flag=='4rcl':
        x15 , x16 , y15, y16 = 0.625, 0.875,  0.25, 0.38
        _ , _, ax = rect_plot(x15,x16,y15,y16,ax,'sink')
    elif flag == '1_rect':
        return_ax= True 
        circle1=plt.Circle((.8, .5), 0.05, color='red', fill=False)
        ax.add_artist(circle1)
    else:
        print('Flag not defined!')
    if return_ax == True:
        return ax


def rect_plot(x1,x2,y1,y2,ax,flag='source'):
    G_source = nx.path_graph(5)
    pos_source = {0: (x1, y1),
                1: (x1, y2),
                2: (x2,y2),
                3: (x2,y1),
                4: (x1,y1)}
    if flag == 'source':
        color='g'
    elif flag == 'sink':
        color='r'
    nx.draw_networkx(G_source, pos_source, node_size=5,width=5,edge_color=color,node_color=color, with_labels = False, ax = ax)
    return G_source, pos_source, ax

test='ye'
if test=='yes':
    folder_name = '1_b12_25dv_spf_11_center_circle/'
    ndiv=25
    flag_source='center'
    flag_sink='circle'
    source_sink_generator('./runs/'+folder_name, ndiv, flag_source, flag_sink)
    source_sink_preprocess('./runs/'+folder_name)
