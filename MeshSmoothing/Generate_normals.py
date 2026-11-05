import pymeshlab
import os

def compute_and_save_normals(input_dir, output_dir):
    """
    Generate and save normals for all mesh files in the input directory.
    Args:
        input_dir (str): Path to the directory containing the input meshes.
        output_dir (str): Path to the directory where the output meshes with normals will be saved.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # List all mesh files in the input directory
    mesh_files = [f for f in os.listdir(input_dir) if f.endswith(('.ply', '.obj', '.stl'))]

    for mesh_file in mesh_files:
        try:
            # Initialize a PyMeshLab mesh set
            ms = pymeshlab.MeshSet()

            # Load the mesh
            input_path = os.path.join(input_dir, mesh_file)
            ms.load_new_mesh(input_path)

            # Compute normals
            ms.compute_normal_for_point_clouds()
            ms.compute_normal_per_face()

            # Save the mesh with normals
            output_path = os.path.join(output_dir, mesh_file)
            ms.save_current_mesh(output_path, save_vertex_normal=True)

            print(f"Processed and saved: {mesh_file}")
        except Exception as e:
            print(f"Error processing {mesh_file}: {e}")

# Example usage
input_directory = 'folder_with_smoothed_meshes'
output_directory = 'output_location'
compute_and_save_normals(input_directory, output_directory)
