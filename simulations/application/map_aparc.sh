#!/bin/bash

datadir=../data/application
meshdir=$datadir/mesh
ID=fsl_mni152_2mm
meshfile=$meshdir/$ID.brain.32.mesh
aparc=$meshdir/aparc+aseg.mgz

source activate fenicsproject
python3 python/convert_to_dolfin_mesh.py --meshfile $meshfile --hdf5file $meshdir/$ID.h5 --xdmfdir $meshdir/$ID
python3 python/add_parcellations.py --in_parc $aparc --in_hdf5 $meshdir/$ID.h5 --out_hdf5 $meshdir/$ID-mapped.h5 
python3 python/map_parcellation.py --in_aparc $aparc --in_hdf5 $meshdir/$ID.h5 --out_xdmf $meshdir/$ID-parc.xdmf
conda deactivate
