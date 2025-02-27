# Script made editing script of:
# Agnese Marcato - agnese.marcato@polito.it
# script for the sedimentation of a cloud of spheres
# used for the dataset creation in https://pubs.acs.org/doi/10.1021/acs.iecr.1c04760 
# tutorial example: https://gitlab.com/yade-dev/trunk/blob/master/doc/sphinx/tutorial/02-gravity-deposition.py
# For more infos check: https://github.com/mulmopro/YadePacking

# Python's libraries to import

### Create a folder for each simulation
def create_folder(dir_name):
	folder = os.getcwd()+'/'+dir_name
	if os.path.exists(folder)==False: 
		os.mkdir(folder)

### Particle size distribution
def particle_size_distribution(N_particles, x_axis_distribution, y_axis_pdf, dir_name, want_plot_pdf=True, want_plot_cdf=True):
	#N_particles = 4284	   # number of particles
	if x_axis_distribution.shape!=y_axis_pdf.shape:
		print('x and y of the distribution must be of the same length :(')
	y_axis_pdf = np.round((y_axis_pdf*N_particles),0).astype(int) # (integer) number of particles for every class
	if want_plot_pdf:
		plot_pdf(x_axis_distribution,y_axis_pdf,dir_name)
	if want_plot_cdf:
		plot_cdf(x_axis_distribution,y_axis_pdf,dir_name)
	return y_axis_pdf

# plot of the resulting distribution (PDF)
def plot_pdf(x_axis_distribution, y_axis_pdf, dir_name):
	plt.figure(figsize=[6,6])
	plt.bar(x_axis_distribution,y_axis_pdf) # new distribution
	plt.xlabel(r'Diameter [$\mu m$] ')
	plt.ylabel('Count')
	plt.savefig(dir_name+'/psd.png')

# plot of the resulting distribution (CDF)
def plot_cdf(x_axis_distribution, y_axis_pdf, dir_name):
	y_axis_cdf=np.zeros(len(y_axis_pdf))
	for i in range(0,len(y_axis_pdf)):
		y_axis_cdf[i] = sum(y_axis_pdf[:(i+1)]) 
	y_axis_cdf=y_axis_cdf/(max(y_axis_cdf))
	plt.figure(figsize=[6,6])
	plt.plot(x_axis_distribution,y_axis_cdf)
	plt.xlabel(r'Diameter [$\mu m$] ')
	plt.ylabel('Count')
	plt.savefig(dir_name+'/cdf.png')

def find_location_in_mesh(Cx,Cy,Cz,radius,x_min,x_max,y_min,y_max,z_min,z_max,dir_name, output_filename):
	''' Find point inside box and outside spheres: for OpenFOAM simulations '''
	size = x_max - x_min
	
	guess = np.random.rand(3)*size
	guess[0] += x_min
	guess[1] += y_min
	guess[2] += z_min
	print(guess)
	# number of trials
	iterates = 1000

	for i in range(0,iterates):
		# distance between centre of the spheres and guess point
		distance = np.sqrt(np.power(Cx-guess[0],2)+np.power(Cy-guess[1],2)+np.power(Cz-guess[2],2))
		#comparison between distance and radius for every sphere
		compare = np.greater(radius, distance*0.9)
		# if radius is never greater than distance we found the point
		result = np.any(compare)
		if result == False:
			fileOut=open(dir_name+"/"+output_filename+".txt","a")
			print('Point location in mesh: (%r,%r,%r)' %(guess[0], guess[1], guess[2]))
			fileOut.write('// Point location in mesh: (%r,%r,%r)\n' %(guess[0], guess[1], guess[2]))
			print('Numero di iterate: %r' %i)
			fileOut.close()		   
			break
		else:
			guess = np.random.rand(3)*size
			guess[0] += x_min
			guess[1] += y_min
			guess[2] += z_min
	return guess

def VTK_creation(Cx,Cy,Cz,radius,dir_name,filename='packing'):
	# Creation of VTK file 
	VTK = open(dir_name+'/'+filename+'.vtk','w')
	VTK.write('# vtk DataFile Version 3.0.\n')
	VTK.write('comment\n')
	VTK.write('ASCII\n\n')
	VTK.write('DATASET POLYDATA\n')
	VTK.write('POINTS %i double\n' %len(Cx))
	for i in range(0,len(Cx)):
		VTK.write('%f %f %f\n' %(Cx[i],Cy[i],Cz[i]))
	VTK.write('\nPOINT_DATA %i\n' %len(Cx))
	VTK.write('SCALARS radius double 1\n')
	VTK.write('LOOKUP_TABLE default\n')
	for i in range(0,len(Cx)):
		VTK.write('%f\n' %(radius[i]))
	VTK.close()

def VTK_blender(Cx,Cy,Cz,radius,dir_name,filename='blender', dim=0):
	# Creation of VTK file 
	VTK = open(dir_name+'/'+filename+'.vtk','w')
	for i in range(0,len(Cx)):
		if(abs(Cx[i])<dim/2 and abs(Cy[i])<dim/2 and Cz[i]>0): 
			VTK.write('%f %f %f %f\n' %(Cx[i],Cy[i],Cz[i],radius[i]))
	VTK.close()

def plot_hist(radius, dir_name, filename='psd'):
	# we check the psd in the bulk cube
	plt.figure(figsize=[6,6])
	plt.hist(radius,bins=int(m.log(len(radius),2)+1)) # Sturge's formula https://en.wikipedia.org/wiki/Histogram
	plt.savefig(dir_name+'/'+filename+'.png')


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

def meanSize(x_axis_distribution, numFrac):
    d4 = 0
    d3 = 0
    for (diam, ni) in zip(x_axis_distribution, numFrac):
        d3 += ni * (diam**3)
        d4 += ni * (diam**4)

    meanSize = d4/d3
    return meanSize

def maxSize(x_axis_distribution, y_axis_pdf):
    mask = np.ones(np.size(y_axis_pdf), dtype=bool)
    for i in range(np.size(y_axis_pdf)):
        if y_axis_pdf[i] == 0:
            mask[i] = False
    x_axis_distribution_clean = x_axis_distribution[mask,...]
    maxSize = x_axis_distribution_clean[-1]
    return maxSize
   
from yade import pack, plot, export
import math as m
import numpy as np
import matplotlib.pyplot as plt
#from functions_script import *
#from sizeCalculation import *
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

Raw_x_axis_distribution = np.array(dati_setup.get("Raw_x_axis_distribution", []))
Raw_numFrac = np.array(dati_setup.get("Raw_numFrac", []))

print("scale:", scale)
print("L:", L)
print("Lz:", Lz)
print("N_particles:", N_particles)
print("Raw_x_axis_distribution:", Raw_x_axis_distribution[0])
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
