import os
import pymeshlab

# Directories
input_dir = "folder_with_meshes"
output_dir = "output_folder"
os.makedirs(output_dir, exist_ok=True)

for filename in os.listdir(input_dir):
    if filename.endswith(".ply"):
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename.replace(".ply", "_smoothed.ply"))
        
        # Mesh processing
        ms = pymeshlab.MeshSet()
        ms.load_new_mesh(input_path)
        ms.apply_coord_hc_laplacian_smoothing()
        ms.apply_coord_hc_laplacian_smoothing()
        # Sampling
        ms.generate_sampling_poisson_disk(samplenum=20000, montecarlorate=20) #20000,20
        
        # Surface reconstruction
        ms.generate_surface_reconstruction_screened_poisson(depth=8, scale=1.1) #8,1.1
        
        ms.meshing_close_holes(maxholesize=1000)
        
        # Smoothing
        
        ms.apply_coord_hc_laplacian_smoothing()
        ms.apply_coord_hc_laplacian_smoothing()
        ms.apply_coord_hc_laplacian_smoothing()

        ms.meshing_surface_subdivision_loop(
            iterations=10,
            selected=False
        )

        ms.meshing_close_holes(
           maxholesize=1)


        # Save the result
        ms.save_current_mesh(output_path)
        print(f"Processed: {output_path}")
