import ovito
import ovito.pipeline
import scipy
import ase
import math
import scipy
import numpy as np
from ovito.io import import_file, export_file
from ovito.modifiers import FreezePropertyModifier ,CombineDatasetsModifier, ClusterAnalysisModifier, DeleteSelectedModifier, ExpressionSelectionModifier, WignerSeitzAnalysisModifier, InvertSelectionModifier
from ase.visualize import view
from ase import Atoms
from ase.cluster import Octahedron
from ovito.io.ase import ovito_to_ase, ase_to_ovito
from ovito.pipeline import Pipeline, StaticSource

def find_atom_cluster(input_file: str, frame: int = -1, cluster_num: int = 1, max_size: int = 1,analysis_method:str = "WS", cutoff: int = 3.2, pot_e: int = -163):
    
    lattice_constant = 5.468
    
    max_leng = max_size * lattice_constant/2

    pipeline = import_file(input_file)
    match analysis_method:
        case "WS":
            pipeline.modifiers.append(WignerSeitzAnalysisModifier(output_displaced=True, reference=pipeline.source))
            pipeline.modifiers.append(ExpressionSelectionModifier(expression='Occupancy > 1 || Occupancy == 0'))

        case "PE":
            pipeline.modifiers.append(ExpressionSelectionModifier(expression='PotentialEnergy >='+str(pot_e)))

    pipeline.modifiers.append(InvertSelectionModifier())
    pipeline.modifiers.append(DeleteSelectedModifier())
    pipeline.modifiers.append(ClusterAnalysisModifier(cluster_coloring=True, cutoff=cutoff, only_selected=False, sort_by_size=True, compute_com=True))
    res_data = pipeline.compute(frame)
    
    com_arr = res_data.tables['clusters']
    com_arr = np.array(com_arr['Center of Mass'])
    com = com_arr[cluster_num -1]


    part = res_data.particles.cluster[...]
    pos_original = res_data.particles.positions_
    filter_arr = []
    for p in part:  # isolate positions for 1 cluster
        if p == cluster_num:
            filter_arr.append(True)
        else:
            filter_arr.append(False)
    pos = pos_original[filter_arr]
    
    new_pos = np.array(pos)
    size = len(new_pos)
    for i in range(size):
        new_pos[i] = new_pos[i]-com
    
    length_x = max(abs(new_pos[:,0]))
    length_y = max(abs(new_pos[:,1]))
    length_z = max(abs(new_pos[:,2]))
    cube_length = max([length_x,length_y,length_z])  #cube side length /2 minimum cube length, smallest cube can be
    cube_radius = math.sqrt(3) * cube_length
    
    
    #find the D
    inv_filter_arr = []
    for b in filter_arr: inv_filter_arr.append( not b)
    not_in_cluster =  pos_original[inv_filter_arr]
    pos_around_com = np.array(not_in_cluster)
    for j in range(len(not_in_cluster)):
        pos_around_com[j] =pos_around_com[j] -com

    lengths = np.zeros(len(not_in_cluster)) #to particles not in cluster
    for j in range(len(not_in_cluster)):
        lengths[j] = math.sqrt((pos_around_com[j])[0]**2 + (pos_around_com[j])[1]**2+(pos_around_com[j])[2]**2 )
    
    lengths = np.sort(lengths)
    count = 0
    while count < len(lengths)-1 and lengths[count] < cube_radius:
        count += 1
    d = 0
    if lengths.size == 0: # for the specific case when there is only 1 cluster
        d = max_leng * math.sqrt(3)
    else:
        d = lengths[count]

    optimal_cube_length = d/math.sqrt(3) #biggest the cube can be
    used_cube_length = min([optimal_cube_length,max_leng])
    
    #the value for the side of the cube, has to be in between optimal_cub_length and cube_length
    pipeline2 = import_file(input_file)
    expr = 'Position.X > '+str(com[0]-used_cube_length)+' && Position.X < '+str(com[0]+used_cube_length)+' && Position.Y > '+str(com[1]-used_cube_length)+' && Position.Y <'+str(com[1]+used_cube_length)+' && Position.Z >'+str(com[2]-used_cube_length)+' && Position.Z <'+str(com[2]+used_cube_length)
    
    
    pipeline2.modifiers.append(ExpressionSelectionModifier(expression=expr))
    pipeline2.modifiers.append(InvertSelectionModifier())
    pipeline2.modifiers.append(DeleteSelectedModifier())
    #if used_cube_length < max_leng then we combine cells:
    
    if used_cube_length !=max_leng:
            #only for combining cells
        pipeline3 = import_file(input_file)
        pipeline3.modifiers.append(ExpressionSelectionModifier(expression=expr))
        pipeline3.modifiers.append(DeleteSelectedModifier())
        new_cell = pipeline3.compute(0)
        export_file(new_cell, "new_cell_for_bigcluster", 'xyz',columns =
                ["Particle Identifier", "Particle Type","velocities", "Force", "Potential Energy",  "Position.X", "Position.Y", "Position.Z"])

        modifer = CombineDatasetsModifier()
        modifer.source.load("new_cell_for_bigcluster")
        pipeline2.modifiers.append(modifer)
        
        #for deleting overlaps due to vibrations
        ws = WignerSeitzAnalysisModifier(output_displaced=True, reference_frame=0)
        ws.reference = ovito.pipeline.FileSource()
        ws.reference.load(input_file)
        pipeline2.modifiers.append(ws)
        pipeline2.modifiers.append(ExpressionSelectionModifier(expression='Occupancy > 1 && ParticleIdentifier != 0'))
        pipeline2.modifiers.append(DeleteSelectedModifier())

        expr_max = 'Position.X > '+str(com[0]-max_leng)+' && Position.X < '+str(com[0]+max_leng)+' && Position.Y > '+str(com[1]-max_leng)+' && Position.Y <'+str(com[1]+max_leng)+' && Position.Z >'+str(com[2]-max_leng)+' && Position.Z <'+str(com[2]+max_leng)
    
        pipeline2.modifiers.append(ExpressionSelectionModifier(expression=expr_max))
        pipeline2.modifiers.append(InvertSelectionModifier())
        pipeline2.modifiers.append(DeleteSelectedModifier())

        
    cluster = pipeline2.compute(frame=frame)
    

    
    
    return cluster


def num_clusters(input_file: str, frame: int):

    pipe = import_file(input_file)

    pipe.modifiers.append(WignerSeitzAnalysisModifier(output_displaced=True, reference=pipe.source))
    pipe.modifiers.append(ExpressionSelectionModifier(expression='Occupancy ==1'))
    pipe.modifiers.append(DeleteSelectedModifier())
    pipe.modifiers.append(ClusterAnalysisModifier(cluster_coloring=True, cutoff=3.2, only_selected=False, sort_by_size=True))

    data = pipe.compute(frame)
    res = data.attributes['ClusterAnalysis.cluster_count']
    return res


#numOfCluster =  #need to figure out cluster_analysis_1 has the exact pipeline to calculate this value, so I can just use that
# for j in range(1, numOfCluster):       # for testing all clusters
"""
for i in range(1,21): # is for testing multiple values

    
    
    cluster_counter = num_clusters("200eV/res-PKA-"+str(i)+"/trajectory_out.xyz", 100)
    for j in range(1,cluster_counter+1):
        test = find_atom_cluster("200eV/res-PKA-"+str(i)+"/trajectory_out.xyz", frame=100,cluster_num=j, max_size=3)
        print(test)
        export_file(test, "200eV test files/200eV_test3_cluster_"+str(j)+"_extract_sample"+str(i), 'xyz',columns =
            ["Particle Identifier", "Particle Type","velocities", "Force", "Potential Energy",  "Position.X", "Position.Y", "Position.Z"])
"""
test = find_atom_cluster(input_file = "1kev/res-PKA-1/trajectory_out.xyz", frame=100, cluster_num = 1,analysis_method="PE", max_size =3,cutoff = 4, pot_e = -162.75)

print(test)
export_file(test, "1keV test files/1keV_test5_cluster_1_extract_sample1", 'xyz',columns =
            ["Particle Identifier", "Particle Type","velocities", "Force", "Potential Energy",  "Position.X", "Position.Y", "Position.Z"])












    

"""
#maybe use it if I need to specify size of cell ,where there is close clusters.
def simple_extraction(input_file: str, size: int , center:[int,int,int]):
    num_atoms = size*size*size*8
    si_lattice_a = 5.468
    pipeline = import_file(input_file)
    data = pipeline.compute(0)

    amountpos = len(data.particles.positions)
    
    vol = data.attributes['volume']
    sideL = math.cbrt(vol)
    amountPerSide = math.cbrt(amountpos)
    leng = sideL/amountPerSide
    leng = si_lattice_a/2

    expr = 'Position.X >= '+str(center[0]-leng*size)+' && Position.X <= '+str(center[0]+leng*size)+' && Position.Y >= '+str(center[1]-leng*size)+' && Position.Y <='+str(center[1]+leng*size)+' && Position.Z >='+str(center[2]-leng*size)+' && Position.Z <='+str(center[2]+leng*size)
    

    pipeline2 = import_file(input_file)
    pipeline2.modifiers.append(ExpressionSelectionModifier(expression=expr))
    pipeline2.modifiers.append(InvertSelectionModifier())
    pipeline2.modifiers.append(DeleteSelectedModifier())
    res = pipeline2.compute(0)

    resu = len(res.particles.positions)
    return (resu,num_atoms,res)

#test = simple_extraction("trj-4.xyz",2)
#print(test)
#export_file(test[2], "extraction_test", 'xyz',columns =
#        ["Particle Identifier", "Particle Type","velocities", "Force", "Potential Energy",  "Position.X", "Position.Y", "Position.Z"])
"""
