READ_ME


Code defect_cluster_isolator.py uses OVITO and ase to extract a singular defect cluster:

variables:
	- input_file: the PATH to the target trajectory file

	- frame: which frame should the cluster be extracted from. Default is last frame.

	- cluster_num: Which cluster should be extracted by number, in order of size, so 1 is biggest and so on.

	- max_size: The amount of unit cells that should be extracted around the defect cluster. so 1 would be 8 atoms, 2 would be 64 atoms. calculated with x^3*8

	- analysis_method: How the defect clusters should be identified. "WS" will use wigner seits defect analysis, and "PE" will use potential energy with the pot_e variable cutoff

	- cutoff: the cutoff length for the cluster analysis modifier.

	- pot_e: is only used if analysis_method is on "PE", and defines the potential energy cutoff value.


The code utilizes 4 pipelines from ovito to isolate the desired size box of atoms with the target defect cluster. Then It utilizes ASE to fix the simulation cell to match the desired size.

1. pipeline is a data extraction. first using which ever analysis specified earlier, to identify the defects. Then deleting all non-defective atoms, and then clustering the defects.

The data collected from pipeline 1., is used to identify the positions of the Center of Mass (COM) of the target defect cluster and the positions of each defective particle. Then it finds the from the COM to all particles not in the target defect cluster.
The it find the closest particle to the target defect cluster, which is outside the defective cluster. call this distance optimal_cube_length,and take the minimum of optimal_cube_length and max_leng, as the used_cube_length. 
Then it is locked into the grid between the atoms. To create new_used_cube_length. Then the expression for selection the target cube with the defect cluster in, is made from this value.

If the used_cube_length == max_leng then the cube is directly extracted in pipeline 2., if not then pipeline 3. is used to create a perfect cell from used_cube_length up to the target max_leng. By extracting the cube made by new_used_cube_length from last frame and extracting the inverse from last frame, then combining the two files and extracting a cube with dimensions from max_leng.


At last the data frame pipeline 2. is made into an Atoms object from ASE, and the simulation cell is changed to the desired dimensions from max_size. To export more nicely a pipeline 4. from ovito is made. 

The returned result is the isolated cube with the target defect inside a perfect cell of desired size.



Code extraction_of_reference_image.py for Wigner Seits analysis:

To run Wigner Seits analysis on the result from the code defect_cluster_isolator.py, you will need a reference image. this can be made from the exact same code, just change the compute from pipeline 2. to compute frame 0.

change: 	cluster = pipeline2.compute(frame=frame_calc) --> cluster = pipeline2.compute(frame=0)

then you will get a reference frame for the result gathered in defect_cluster_isolator.py, which is what this code file is.