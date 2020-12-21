import sys
#-----------------------------
path2core = '../img2net-core/'
if path2core not in sys.path:
    sys.path.append(path2core)
import img2net


image_path = 'input/angolan_river_crop.png'

img2net.image2net(image_path, 
	N_runs=3, 
	t2 = .25, 
	t3 = .6, 
	new_size = 300, 
	reversed_colors =True, 
	weighting_method_simplification = 'ER', 
	ds_interpolation = 'area')