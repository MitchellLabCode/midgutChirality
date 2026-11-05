import pymeshlab
import os

directory_path= 'smoothed_meshes_location' #Directory with your meshes (.ply files)



#Iterates through all meshes and checks if genus is above 0, in which case it prints out which mesh is the issue and what the genus is



genus_check_list=[]
for filename in os.listdir(directory_path):
    file_path = os.path.join(directory_path, filename)
    if "mesh" in filename:
        ms = pymeshlab.MeshSet()
        ms.load_new_mesh(file_path)
        ms.get_topological_measures()
        topology = ms.get_topological_measures()
        if "genus" in topology and topology["genus"] > 0:
        #print(f"{filename} has genus {topology['genus']}")
            genus_check_list.append(filename)

    genus_check_list.sort()
print(genus_check_list)
print(f"Number of files with genus greater than 0: {len(genus_check_list)}")
