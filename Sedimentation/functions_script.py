import os
import math as m
import numpy as np
import matplotlib.pyplot as plt


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

