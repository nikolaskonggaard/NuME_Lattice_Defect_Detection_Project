import ovito
import ovito.pipeline

import math
import scipy
import numpy as np
from ovito.io import import_file, export_file
from ovito.modifiers import CombineDatasetsModifier, ClusterAnalysisModifier, DeleteSelectedModifier, ExpressionSelectionModifier, WignerSeitzAnalysisModifier, InvertSelectionModifier
from ase import Atoms
from ovito.io.ase import ovito_to_ase, ase_to_ovito
from ovito.pipeline import Pipeline, StaticSource


#inputs are the file name, frame number default last frame, cluster number from biggest to smallest: 1=biggest.
#analysis method: WS= wignerseits, PE = Potential energy
#cutoff distance for cluster analysis, and pot_e which is the potential energy cutoff, only necessary if analysis_method="PE"
def find_atom_cluster(input_file: str, frame: int = -1, cluster_num: int = 1, max_size: int = 1,analysis_method:str = "WS", cutoff: int = 3.2, pot_e: int = -163):
    
    lattice_constant = 5.468
    
    
    max_leng = max_size * lattice_constant/2
    

    pipeline = import_file(input_file)
    start_data = pipeline.source.compute()

    #create pipeline to find data of target cluster, depending on which analysis method is used.
    lattice_unit = start_data.cell[0][0]/((len(start_data.particles_.positions) /8)**(1/3) * 2)
    match analysis_method:
        case "WS":
            pipeline.modifiers.append(WignerSeitzAnalysisModifier(output_displaced=True, reference=pipeline.source))
            pipeline.modifiers.append(ExpressionSelectionModifier(expression='Occupancy > 1 || Occupancy == 0'))

        case "PE":
            pipeline.modifiers.append(ExpressionSelectionModifier(expression='PotentialEnergy >='+str(pot_e)))

    pipeline.modifiers.append(InvertSelectionModifier())
    pipeline.modifiers.append(DeleteSelectedModifier())
    pipeline.modifiers.append(ClusterAnalysisModifier(cluster_coloring=True, cutoff=cutoff, only_selected=False, sort_by_size=True, compute_com=True))
    
    #choose which frame will be used
    #negative values and values outside range of frame in trajectory file will give the last frame
    frame_calc = frame
    if frame < 0 | frame > pipeline.source.num_frames:
        frame_calc = pipeline.source.num_frames
        
    #data to use in calculation
    res_data = pipeline.compute(frame_calc)
    
    #Center of Mass
    com_arr = res_data.tables['clusters']
    com_arr = np.array(com_arr['Center of Mass'])
    com = com_arr[cluster_num -1]

    #Here we find the positions of all the particles in our target cluster
    #by making an array of boolean and filter all original atom positions with regards to the origin
    part = res_data.particles.cluster[...]
    pos_original = res_data.particles.positions_
    filter_arr = []
    for p in part:  # isolate positions for 1 cluster
        if p == cluster_num:
            filter_arr.append(True)
        else:
            filter_arr.append(False)
    pos = pos_original[filter_arr]
    
    #find the positions of the target clusters atoms with regards to the center of mass
    new_pos = np.array(pos)
    size = len(new_pos)
    for i in range(size):
        new_pos[i] = new_pos[i]-com
    
    #find the length of the atom in the target cluster that is furthest away from the center of mass
    length_x = max(abs(new_pos[:,0]))
    length_y = max(abs(new_pos[:,1]))
    length_z = max(abs(new_pos[:,2]))
    cube_length = max([length_x,length_y,length_z])  #cube side length /2 minimum cube length, smallest cube can be
    cube_radius = math.sqrt(3) * cube_length         #makes sure that the cube have all atos inside it
    
    
    #find the Distance to the closest defect clusters that are not the target
    inv_filter_arr = []
    for b in filter_arr: inv_filter_arr.append( not b)
    not_in_cluster =  pos_original[inv_filter_arr]
    pos_around_com = np.array(not_in_cluster)
    for j in range(len(not_in_cluster)):
        pos_around_com[j] =pos_around_com[j] -com

    lengths = np.zeros(len(not_in_cluster)) #to particles not in cluster
    for j in range(len(not_in_cluster)):
        lengths[j] = math.sqrt((pos_around_com[j])[0]**2 + (pos_around_com[j])[1]**2+(pos_around_com[j])[2]**2 )
    
    #filters out the defects that are closer to the center of mass of the target defect cluster than
    #atoms in the target cluster itself. 
    #Therefore finding the nearest defect outside the cluster, to the cluster
    lengths = np.sort(lengths)
    count = 0
    while count < len(lengths)-1 and lengths[count] < cube_radius:
        count += 1
    d = 0
    if lengths.size == 0: # for the specific case when there is only 1 cluster
        d = max_leng * math.sqrt(3)
    else:
        d = lengths[count]

    #Finding the length we should extract from the last frame, which has the whole defect.
    #The minimum between the target length (max_leng) and the closest defect which is not part of the cluster (optimal_cube_length)
    optimal_cube_length = d/math.sqrt(3) #biggest the cube can be
    used_cube_length = min([optimal_cube_length,max_leng])
    
    #Locking the box to the lattice grid, to make sure nice boundaries are created.
    max_used_edges = [com[0]+used_cube_length, com[1]+used_cube_length, com[2]+used_cube_length] #the max of the cubes
    min_used_edges = [com[0]-used_cube_length, com[1]-used_cube_length, com[2]-used_cube_length] # the min of edges
    new_max_used_edges = []
    new_min_used_edges =[]
    for ed in max_used_edges:
        i = 0
        while ed > lattice_unit * (1/4 + 1/2 * (i)):
            i += 1
        new_max_used_edges.append(lattice_unit * (1/4 + 1/2 * i))
    for ed in min_used_edges:
        i = 0
        while ed > lattice_unit * (1/4 + 1/2 * (i)):
            i += 1
        new_min_used_edges.append(lattice_unit * (1/4 + 1/2 * (i-1)))

    #Create a new pipeline for extracting the cluster, using the length we just calculated.
    pipeline2 = import_file(input_file)

    #expression for the cube with the target cluster
    expr = 'Position.X > '+str(new_min_used_edges[0])+' && Position.X < '+str(new_max_used_edges[0])+' && Position.Y > '+str(new_min_used_edges[1])+' && Position.Y <'+str(new_max_used_edges[1])+' && Position.Z >'+str(new_min_used_edges[2])+' && Position.Z <'+str(new_max_used_edges[2])
    
    #Find the target box using the center of mass of target cluster and the intended length (max_leng)
    #Here it is locking it to the grid again for the new bigger box
    max_edges = [com[0]+max_leng, com[1]+max_leng, com[2]+max_leng] #the max of the cubes
    min_edges = [com[0]-max_leng, com[1]-max_leng, com[2]-max_leng] # the min of edges
    new_max_edges = []
    new_min_edges =[]
    for ed in max_edges:
        i = 0
        while ed > lattice_unit * (1/4 + 1/2 * (i)):
            i += 1
        new_max_edges.append(lattice_unit * (1/4 + 1/2 * i))
    for ed in min_edges:
        i = 0
        while ed > lattice_unit * (1/4 + 1/2 * i):
            i += 1
        new_min_edges.append(lattice_unit * (1/4 + 1/2 * (i))) #have i-1 if bigger cube
        
    #expression for the final cube
    new_expr = 'Position.X > '+str(new_min_edges[0])+' && Position.X < '+str(new_max_edges[0])+' && Position.Y > '+str(new_min_edges[1])+' && Position.Y <'+str(new_max_edges[1])+' && Position.Z >'+str(new_min_edges[2])+' && Position.Z <'+str(new_max_edges[2])
    
    #checking if the cube with the target cluster is the intended size or smaller
    #If smaller then we have to isolate it and put it into frame 0, to extract a perfect buffer zone around
    if used_cube_length !=max_leng:

        pipeline2.modifiers.append(ExpressionSelectionModifier(expression=expr))
        pipeline2.modifiers.append(InvertSelectionModifier())
        pipeline2.modifiers.append(DeleteSelectedModifier())
            #only for combining cells, it creates a file for the perfect cell,
            # that has a hole where the cluster will be inserted
        pipeline3 = import_file(input_file)
        pipeline3.modifiers.append(ExpressionSelectionModifier(expression=expr))
        pipeline3.modifiers.append(DeleteSelectedModifier())
        new_cell = pipeline3.compute(0)
        export_file(new_cell, "new_cell_for_bigcluster", 'xyz',columns =
                ["Particle Identifier", "Particle Type","velocities", "Force", "Potential Energy",  "Position.X", "Position.Y", "Position.Z"])

        modifer = CombineDatasetsModifier()
        modifer.source.load("new_cell_for_bigcluster")
        pipeline2.modifiers.append(modifer)
        
        #for deleting overlaps due to vibrations, if some atoms overlap when combining the file
        #this shouldnt happen often (or at all), but it is a safety.
        ws = WignerSeitzAnalysisModifier(output_displaced=True, reference_frame=0)
        ws.reference = ovito.pipeline.FileSource()
        ws.reference.load(input_file)
        pipeline2.modifiers.append(ws)
        pipeline2.modifiers.append(ExpressionSelectionModifier(expression='Occupancy > 1 && ParticleIdentifier != 0'))
        pipeline2.modifiers.append(DeleteSelectedModifier())

        
    #Isolate the intended cell with a singular defect cluster
    pipeline2.modifiers.append(ExpressionSelectionModifier(expression=new_expr))
    pipeline2.modifiers.append(InvertSelectionModifier())
    pipeline2.modifiers.append(DeleteSelectedModifier())


    cluster = pipeline2.compute(frame=frame_calc)


#ase polishing
    #using ase change the simulation cell to have pbc and fit nicely around the cell
    polishing = Atoms(symbols=str(len(cluster.particles.positions))+"Si", positions=cluster.particles.positions_, pbc=True)
    polishing.center(vacuum=0)
    
    extended = max_size * lattice_constant
    mod_cell = [[extended,0,0],[0,extended,0],[0,0,extended]]
    polishing.set_cell(cell=mod_cell)
    
    
    
    #Make it into ovito data for exporting it easilier 
    ase_ovito_data = ase_to_ovito(polishing)
    pipeline4 = Pipeline(source = StaticSource(data=ase_ovito_data))
    res = pipeline4.compute()
    
    return res



test = find_atom_cluster(input_file = "PATH/file.xyz", cluster_num = 1,analysis_method="WS", max_size =5,cutoff = 3.2, pot_e = -162.75)
export_file(test, "PATH/Result.xyz", 'xyz',columns =
           ["Particle Identifier", "Particle Type",  "Position.X", "Position.Y", "Position.Z"]) 
