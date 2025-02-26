cd Sedimentation
singularity run of9_yade22_py310_latest.sif yade-batch --job-threads 10 sedimentation_pack.py
rm blender99.vtk
cp ./packing/packing_total.vtk ./blender99.vtk
blender -b -P Imp.py
rm b1*
rm -r ./__pycache__
rm -r ./packing
rm sedimentation_pack.py.default.log
cp dim_pos_rot.txt ../geometry.txt
