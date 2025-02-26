# Script made editing script of:
# Agnese Marcato - agnese.marcato@polito.it
# script for the sedimentation of a cloud of spheres
# used for the dataset creation in https://pubs.acs.org/doi/10.1021/acs.iecr.1c04760 
# tutorial example: https://gitlab.com/yade-dev/trunk/blob/master/doc/sphinx/tutorial/02-gravity-deposition.py
# For more infos check: https://github.com/mulmopro/YadePacking

# Python's libraries to import

def leggi_setup(file_path):
    dati = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("#"):  # Ignora righe vuote o commentate
                continue
            chiave, valore = map(str.strip, line.split("=", 1))
            try:
                dati[chiave] = ast.literal_eval(valore)
            except (ValueError, SyntaxError):
                dati[chiave] = valore
    return dati

from yade import pack, plot, export
import math as m
import numpy as np
import matplotlib.pyplot as plt
from functions_script import *
from sizeCalculation import *
import inspect, os
import ast
frame = inspect.currentframe()
filePath = inspect.getfile(frame)
currentDir = os.path.realpath(os.path.abspath(os.path.dirname(filePath)))




####### USER INPUTS 
dir_name='packing'    # name of the directory containing all the simulation files
want_plot_cdf = True        # do you want to plot the pdf of the distribution?
want_plot_pdf = True        # do you want to plot the cdf of the distribution? 
output_filename='Outfile'

create_folder(dir_name)

dati_setup = leggi_setup(currentDir+"/../setup.txt")
scale = dati_setup.get("scale", 1) 
L = dati_setup.get("L", 5)
Lz = dati_setup.get("Lz", 10 * L)
N_particles = dati_setup.get("N_particles", 100)

Raw_x_axis_distribution = np.array(dati_setup.get("Raw_x_axis_distribution", [])) / scale
Raw_numFrac = np.array(dati_setup.get("Raw_numFrac", []))

print("scale:", scale)
print("L:", L)
print("Lz:", Lz)
print("N_particles:", N_particles)
print("Raw_x_axis_distribution:", Raw_x_axis_distribution)
print("Raw_numFrac:", Raw_numFrac)


x_axis_distribution=Raw_x_axis_distribution
numFrac = Raw_numFrac/sum(Raw_numFrac)


meanValue = meanSize(x_axis_distribution, numFrac)

porosity = 0.4              # Hypothetical porosity of the packing

y_axis_pdf = particle_size_distribution(N_particles, x_axis_distribution, numFrac, dir_name, want_plot_pdf=True, want_plot_cdf=True) # let's scale everything

# For L and Lz insert an integer number That will be equal to a multiple of the mean size
#L = box_dimension(x_axis_distribution, y_axis_pdf, porosity)                      # Lx=Ly dimension of the simulation box 


width_box = L/10              # If you want to cut a subsampling set the length of the cube, if you don't want to, just neglect this

###################




####### PACKING CREATION
O.bodies.append(geom.facetBox(
    (0, 0,(Lz+100)/2),     # Center of the box
    (L/2,L/2,(Lz+100)/2),     # Length of the box sides, I usually make the box longer you can change 100
    wallMask=31),
)

sp=pack.SpherePack()
# create empty sphere packing
# sphere packing is not equivalent to particles in simulation, it contains only the pure geometry

fileOut = open(dir_name+"/"+output_filename+".txt","w")
fileOut.write('\nMean Size (micro_m): %.3f' % (meanValue))
fileOut.write('\nDimension of the box: %d' % (L))
fileOut.write('\nHeight of the box: %d\n' % (Lz))
fileOut.close()



for i in range(0,len(x_axis_distribution)): # make a cloud for every class of particles 
    sp.makeCloud(
        ((-L/2),(-L/2),0),((L/2),(L/2),Lz),                                                    # Corners of the dispersion
        rMean    = x_axis_distribution[len(x_axis_distribution)-1-i]/2,      # Mean radius
        #rRelFuzz = 0,                                                       # relative fuzz of the radius 
        num      = int(y_axis_pdf[len(x_axis_distribution)-1-i]),            # Number of spheres
        ) 

sp.toSimulation()

# the simulation loop consists in running defined sequence of 'engines'
O.engines=[
    ForceResetter(),
    InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()]),
    InteractionLoop(
        # handle sphere+sphere and facet+sphere collisions
        [Ig2_Sphere_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
        [Ip2_FrictMat_FrictMat_FrictPhys()],
        [Law2_ScGeom_FrictPhys_CundallStrack()]
    ),
    NewtonIntegrator(gravity=(0,0,-9.81),damping=0.4),
    #qt.SnapshotEngine(fileBase='3d-',iterPeriod=200,label='snapshot'),
    # call the checkUnbalanced function (defined below) every 2 seconds
    PyRunner(command='checkUnbalanced()',realPeriod=2),
]

# definition of the timestep
O.dt=.5*PWaveTimeStep()
O.trackEnergy=True

# if the unbalanced forces goes below the desired value, the packing
# is considered stabilized, therefore stop collect data history and stop.
def checkUnbalanced():
    if unbalancedForce()<0.001:
        print("\nSIMULATION FINISHED")
        O.pause()

O.saveTmp()

O.run()               # run forever, until stopped by checkUnbalanced()

utils.waitIfBatch()

######### POSTPROCESSING

fileOut = open(dir_name+"/"+output_filename+".txt","a")
Cx_box = []
Cy_box = []
Cz_box = []
radius_box = []
Cx_total = []
Cy_total = []
Cz_total = []
radius_total = []
print('Start writing coordinates')
# We want to extract a bulk cube from the entire geometry - we define the x, y, z min and max boundaries of the cube
bulk_cube=18
x_lim_min = -bulk_cube
x_lim_max = bulk_cube
y_lim_min = -bulk_cube
y_lim_max = bulk_cube
zMax = np.max([b.state.pos[2]+b.shape.radius for b in O.bodies if isinstance(b.shape,Sphere)])
z_lim_min = zMax/2-bulk_cube
z_lim_max = zMax/2+bulk_cube


for sph in O.bodies: 
    if isinstance(sph.shape,Sphere):
        Cx_total.append(sph.state.pos[0])
        Cy_total.append(sph.state.pos[1])
        Cz_total.append(sph.state.pos[2])
        radius_total.append(sph.shape.radius)
        if ((sph.state.pos[0]+sph.shape.radius>x_lim_min) & (sph.state.pos[0]-sph.shape.radius<x_lim_max) & (sph.state.pos[1]+sph.shape.radius>y_lim_min) & (sph.state.pos[1]-sph.shape.radius<y_lim_max) & (sph.state.pos[2]+sph.shape.radius>z_lim_min) & (sph.state.pos[2]-sph.shape.radius<z_lim_max)): 
            Cx_box.append(sph.state.pos[0])
            Cy_box.append(sph.state.pos[1])
            Cz_box.append(sph.state.pos[2])
            radius_box.append(sph.shape.radius)
    

# fileOut.write('x_min = 0.0\n')
# fileOut.write('x_max = %.1f\n' % L)
# fileOut.write('y_min = 0.0\n')
# fileOut.write('y_max = %.1f\n' % L)
# fileOut.write('z_min = 0.0\n')
# fileOut.write('z_max = %.1f\n' % (np.max(np.asarray(Cz_box)+np.asarray(radius_box))-1))
# fileOut.close()
fileOut.write('x_min=%r\n' %np.min(np.asarray(Cx_box)-np.asarray(radius_box)))
fileOut.write('x_max=%r\n' %np.max(np.asarray(Cx_box)+np.asarray(radius_box)))
fileOut.write('y_min=%r\n' %np.min(np.asarray(Cy_box)-np.asarray(radius_box)))
fileOut.write('y_max=%r\n' %np.max(np.asarray(Cy_box)+np.asarray(radius_box)))
fileOut.write('z_min=%r\n' %np.min(np.asarray(Cz_box)-np.asarray(radius_box)))
fileOut.write('z_max=%r\n' %np.max(np.asarray(Cz_box)+np.asarray(radius_box)))
fileOut.close()

print('Finish write coordinates')

''' Find point inside box and outside spheres: for OpenFOAM simulations '''
x_lim_min = (-width_box)
x_lim_max = (width_box)
y_lim_min = (-width_box)
y_lim_max = (width_box)
z_lim_min = (zMax/2-width_box)
z_lim_max = (zMax/2+width_box)
print(x_lim_min, x_lim_max, y_lim_min, y_lim_max, z_lim_min, z_lim_max)
guess_point=find_location_in_mesh(Cx_total, Cy_total, Cz_total, radius_total, x_lim_min,x_lim_max,y_lim_min,y_lim_max,z_lim_min,z_lim_max, dir_name, output_filename=output_filename)

print('Finish locate mesh')

plot_hist(radius_total, dir_name, filename='total_psd')

print('Start creating VTK file')
VTK_blender(Cx_total,Cy_total,Cz_total,radius_total,dir_name,filename='packing_total', dim=L)

print('FINISHED\n')
