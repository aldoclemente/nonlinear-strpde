import argparse
import numpy
import nibabel
from nibabel.affines import apply_affine
from dolfin import *

def map_parcellation_to_mesh(parcfile, meshfile, outfile):
    # Load image from the parcellation file,
    # and extract the data it contains
    image  = nibabel.load(parcfile)
    data = image.get_fdata() 

    # Examine the dimensions of the image and
    # examine the tag for the voxel located at 
    # data position 100, 100, 100
    print(data.shape)
    print(data[100, 100, 100])
    # Import brain mesh
    mesh = Mesh()
    
    hdf = HDF5File(mesh.mpi_comm(), meshfile, "r")
    hdf.read(mesh, "/mesh", False)  
    #with XDMFFile(meshfile) as file:
    #    file.read(mesh)
    print(mesh.num_cells())
    
    # Define mesh-based region representation 
    n = mesh.topology().dim()
    regions = MeshFunction("size_t", mesh, n, 0)
    print(regions[0])
    print(regions.array())
    
    # Find the transformation f from T1 voxel space
    # to RAS space and take its inverse to get the
    # map from RAS to voxel space
    vox2ras = image.header.get_vox2ras_tkr()
    ras2vox = numpy.linalg.inv(vox2ras)

    print("Iterating over all cells...")
    for cell in cells(mesh):
        c = cell.index()

        # Extract RAS coordinates of cell midpoint
        xyz = cell.midpoint()[:]

        # Convert to voxel space
        ijk = apply_affine(ras2vox, xyz)

        # Round off to nearest integers to find voxel indices
        i, j, k = numpy.rint(ijk).astype("int")  
        
        # Insert image data into the mesh function:
        regions.array()[c] = int(data[i, j, k])

    # Store regions in XDMF
    xdmf = XDMFFile(mesh.mpi_comm(), outfile)
    xdmf.write(regions)
    xdmf.close()
    
    # ---- Extract arrays ----
    mesh.init()
    nodes = mesh.coordinates().copy()                      # (N,3)
    cells_arr = mesh.cells().copy().astype(numpy.int32)       # (NC,k)

    mesh.init(n-1, n)
    boundary_vertex_ids = []
    for facet in facets(mesh):
        incident_cells = facet.entities(n)
        if len(incident_cells) == 1:                       # exterior facet
            boundary_vertex_ids.extend(facet.entities(0))
    boundary_nodes = numpy.unique(numpy.array(boundary_vertex_ids, dtype=numpy.int32))

    regions_vec = numpy.asarray(regions.array(), dtype=numpy.int32)  # (NC,)

    # ---- Save to text ----
    numpy.savetxt("../data/application/mesh/nodes.txt", nodes, fmt="%.6f")
    numpy.savetxt("../data/application/mesh/cells.txt", cells_arr, fmt="%d")
    numpy.savetxt("../data/application/mesh/boundary.txt", boundary_nodes, fmt="%d")
    numpy.savetxt("../data/application/mesh/regions.txt", regions_vec, fmt="%d")

if __name__ =='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_aparc', type=str)
    parser.add_argument('--in_hdf5', type=str)       
    parser.add_argument('--out_xdmf', type=str)
    Z = parser.parse_args() 
    
    map_parcellation_to_mesh(Z.in_aparc, Z.in_hdf5, Z.out_xdmf)
    
