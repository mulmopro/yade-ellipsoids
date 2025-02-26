def clear():
	bpy.ops.object.select_all(action='SELECT')
	bpy.ops.object.delete()

def deselect():
	bpy.ops.object.select_all(action = 'DESELECT')

def import_stl(name, directory="./", filetype=".stl"):
	bpy.ops.wm.stl_import(filepath = directory + name + filetype)

def change_name(name, newname):
	bpy.context.scene.objects[name].select_set(True)
	bpy.context.view_layer.objects.active = bpy.data.objects[name]
	bpy.context.object.name = newname
	bpy.context.scene.objects[newname].select_set(False)

def position(name, X=None, Y=None, Z=None, dx=0.0, dy=0.0, dz=0.0):
	bpy.context.scene.objects[name].select_set(True)
	bpy.context.view_layer.objects.active = bpy.data.objects[name]
	if X==None:
		X = bpy.context.object.location[0]
	if Y==None:
		Y = bpy.context.object.location[1]
	if Z==None:
		Z = bpy.context.object.location[2]
	bpy.context.object.location = [X+dx, Y+dy, Z+dz]
	bpy.context.scene.objects[name].select_set(False)

def dimension(name, X=None, Y=None, Z=None):
	bpy.context.scene.objects[name].select_set(True)
	bpy.context.view_layer.objects.active = bpy.data.objects[name]
	[A, B, C] = bpy.context.object.dimensions
	if X!=None:
		A=X
	if Y!=None:
		B=Y
	if Z!=None:
		C=Z
	bpy.context.object.dimensions = [A, B, C]
	bpy.context.scene.objects[name].select_set(False)

def rotation(name, X=None, Y=None, Z=None):
	bpy.context.scene.objects[name].select_set(True)
	bpy.context.view_layer.objects.active = bpy.data.objects[name]
	[A, B, C]=bpy.context.object.rotation_euler
	if X!=None:
		bpy.context.object.rotation_euler[0] = X
	if Y!=None:
		bpy.context.object.rotation_euler[1] = Y
	if Z!=None:
		bpy.context.object.rotation_euler[2] = Z
	bpy.context.scene.objects[name].select_set(False)
	return A,B,C
	

def boolean(name1, name2, operation, hole=False):
	bpy.context.scene.objects[name1].select_set(True)
	bpy.context.view_layer.objects.active = bpy.data.objects[name1]
	bpy.ops.object.modifier_add(type = 'BOOLEAN')
	bpy.context.object.modifiers["Boolean"].operation = operation
	bpy.context.object.modifiers["Boolean"].object = bpy.data.objects[name2]
	bpy.context.object.modifiers["Boolean"].use_hole_tolerant = hole
	bpy.ops.object.modifier_apply(modifier = "Boolean")
	bpy.context.scene.objects[name1].select_set(False)

def distance(name, Xc=None, Yc=None, Zc=None, mindist=0.0):
	bpy.context.scene.objects[name].select_set(True)
	bpy.context.view_layer.objects.active = bpy.data.objects[name]
	X = bpy.context.object.location[0]
	Y = bpy.context.object.location[1]
	Z = bpy.context.object.location[2]
	bpy.context.scene.objects[name].select_set(False)
	if Xc==None:
		Xc = X
	if Yc==None:
		Yc = Y
	if Zc==None:
		Zc = Z
	R = np.sqrt((Xc-X)**2 + (Yc-Y)**2 + (Zc-Z)**2)
	if R>mindist:
		return True
	else:
		return False

def scale(name, X=None, Y=None, Z=None, perc=1):
	bpy.context.scene.objects[name].select_set(True)
	bpy.context.view_layer.objects.active = bpy.data.objects[name]
	[A, B, C] = bpy.context.object.scale
	if X!=None:
		A=X
	if Y!=None:
		B=Y
	if Z!=None:
		C=Z
	bpy.context.object.scale = [A*perc, B*perc, C*perc]
	bpy.context.scene.objects[name].select_set(False)

def radial(name, plane="XY"):
	bpy.context.scene.objects[name].select_set(True)
	bpy.context.view_layer.objects.active = bpy.data.objects[name]
	X = bpy.context.object.location[0]
	Y = bpy.context.object.location[1]
	Z = bpy.context.object.location[2]
	bpy.context.scene.objects[name].select_set(False)
	if plane=="XY":
		return np.sqrt(X**2+Y**2)
	if plane=="XZ":
		return np.sqrt(X**2+Z**2)
	if plane=="YZ":
		return np.sqrt(Z**2+Y**2)

def export_stl(filename="output", directory="./", filetype=".stl", ascii_format=True, export_selected_objects=True):
	bpy.ops.object.select_all(action='SELECT')
	bpy.ops.wm.stl_export(filepath=directory+filename+filetype, ascii_format=ascii_format, export_selected_objects=export_selected_objects)

def volume():
	bpy.ops.object.select_all(action='SELECT')
	vol = np.zeros(np.size(bpy.context.selected_objects))
	for i, ob in enumerate(bpy.context.selected_objects):
		matrix = ob.matrix_world
		me = ob.to_mesh()
		me.transform(matrix)
		bm = bmesh.new()
		bm.from_mesh(me)
		vol[i] = bm.calc_volume(signed=True)
		bm.free()

	bpy.ops.object.select_all(action='DESELECT')

	vol = np.round(vol, 4)

	return vol

def check_particles(zmax, numPart, flag):
	obj = bpy.data.objects
	n = 0
	for o in range(np.size(obj)):
		if  obj[o].matrix_world.translation[2] > zmax:
			n += 1
	if n <= (numPart*0.02):
		flag = False

	return flag

def delete(name):
	bpy.context.scene.objects[name].select_set(True)
	bpy.ops.object.delete()
	
def save_state(directory, name):
	bpy.ops.wm.save_as_mainfile(filepath=directory+name)

def roughness(name, rgh=0.1, smt=1):
	bpy.context.scene.objects[name].select_set(True)
	bpy.context.view_layer.objects.active = bpy.data.objects[name]
	bpy.ops.object.editmode_toggle()
	bpy.ops.transform.vertex_random(offset=rgh)
	bpy.ops.mesh.subdivide()
	bpy.ops.mesh.subdivide(smoothness=smt)
	bpy.ops.object.editmode_toggle()
	bpy.context.scene.objects[name].select_set(False)

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

import bpy
import numpy as np
import math, time
import random
import bmesh
import inspect, os
import ast

frame = inspect.currentframe()
filePath = inspect.getfile(frame)
currentDir = os.path.realpath(os.path.abspath(os.path.dirname(filePath)))

startTime = time.time()	

dati_setup = leggi_setup(currentDir+"/../setup.txt")
friction = dati_setup.get("friction", 0.1) 
step_frame = dati_setup.get("step_frame", 100)
max_frame = dati_setup.get("max_frame", 100)
ar = dati_setup.get("Aspect_ratio", 1)
thikness = dati_setup.get("thikness", 3/5)
spheres = dati_setup.get("spheres", True)

#friction = 0.1
restitution = 0
#step_frame = 500
#max_frame = 500
#ar = 0.514
#thikness=3/5
plane="Plane"
rough=False		
rescale=False
square_box=True
#spheres=False



Cyl_Diam=9
f=open("blender99.vtk", 'r')
lines = f.readlines()
f.close
num = np.size(lines)
data = np.zeros((num, 4))
for i, line in enumerate(lines):
	l = line.split()
	for j in range(np.size(l)):
		data[i, j] = float(l[j])

minArr = np.zeros(3)
maxArr = np.zeros(3)
diff = np.zeros(3)
scaleFactor = np.zeros(3)
for i in range(3):
	minArr[i] = np.amin(data[:,i]-data[:,3])
	maxArr[i] = np.amax(data[:,i]+data[:,3])
	diff[i] = maxArr[i] - minArr[i]
	if minArr[i] > 0:
		scaleFactor[i] = 0.99
	else:
		scaleFactor[i] = 1.01
diff = diff.astype(int)
dim = max(diff[0], diff[1]) + 1
print("Lato della box = "+str(dim)+"\n")
dimZ = diff[2] + 1

clear()
out=open("./dim_pos_rot.txt",'w')
DirSTL="./Particle/"
for i in range(num):
	x=data[i, 0]
	y=data[i, 1]
	z=data[i, 2]
	r=data[i, 3]
	obj="p"+str(i)
	import_stl(name="sphere", directory=DirSTL)
	change_name("sphere", obj)
	position(obj, X=(x), Y=(y), Z=(z))
	A=r*2
	if(spheres):
		B=A
		C=A
	else:
		B=A*ar
		C=thikness
		
	dimension(obj, X=A, Y=B, Z=C)
	out.write("D) "+str(A/2*1.02)+" "+str(B/2*1.02)+" "+str(C/2*1.02)+"\n")
	if rough:
		roughness(obj)
f.close()

vol = volume()
m0 = 1/(np.min(vol))
obj = bpy.data.objects
for o in range(np.size(obj)):
	obj[o].select_set(True)
	bpy.ops.rigidbody.objects_add()
	obj[o].rigid_body.type='ACTIVE'
	obj[o].rigid_body.mass = m0*vol[o]
	obj[o].rigid_body.friction = friction
	obj[o].rigid_body.restitution = restitution
	obj[o].rigid_body.collision_shape = 'CONVEX_HULL'
	obj[o].rigid_body.use_deactivation = True
	obj[o].select_set(False)

finishAddObjects = time.time()
print("Time to add objects and define rigid body properties: %.3f" % (finishAddObjects-startTime))

name_list=[]
if(square_box):
	for i in range(4):
		if i<2:
			bpy.ops.mesh.primitive_plane_add(size=dim, enter_editmode=False, align='WORLD', location=(0, ((i*dim)-dim/2), dimZ/2), rotation=(math.pi/2, 0, 0), scale=(1, 1, 1))
			bpy.data.objects[plane].scale[1] = dimZ/dim
			bpy.ops.rigidbody.object_add(type='PASSIVE')
			bpy.context.object.rigid_body.collision_shape = 'MESH'
		else:
			bpy.ops.mesh.primitive_plane_add(size=dim, enter_editmode=False, align='WORLD', location=(((i-2)*dim-dim/2), 0, dimZ/2), rotation=(0, math.pi/2, 0), scale=(1, 1, 1))
			bpy.data.objects[plane].scale[0] = dimZ/dim
			bpy.ops.rigidbody.object_add(type='PASSIVE')
			bpy.context.object.rigid_body.collision_shape = 'MESH'
		name=plane+str(i)
		change_name(plane, name)
		name_list.append(name)
	bpy.ops.mesh.primitive_plane_add(size=dim, enter_editmode=False, align='WORLD', location=(0, 0, 0), rotation=(0, 0, 0), scale=(1, 1, 1))
	bpy.ops.rigidbody.object_add(type='PASSIVE')
	bpy.context.object.rigid_body.collision_shape = 'MESH'
	name=plane+str(4)
	change_name(plane, name)
	name_list.append(name)
else:
	import_stl(name="cylinder")
	dimension("cylinder", X=Cyl_Diam, Y=Cyl_Diam, Z=dimZ)
	bpy.ops.rigidbody.object_add(type='PASSIVE')
	bpy.context.object.rigid_body.collision_shape = 'MESH'
	name_list.append("cylinder")
	
finishAddBox = time.time()
print("Time to add box: %.3f" % (finishAddBox - finishAddObjects))

goNext = True
end_frame = 0
cnt=1
while goNext and (cnt <= max_frame):
	end_frame += step_frame
	bpy.context.scene.frame_end = end_frame
	print('End frame is: '+str(end_frame))
	while (cnt<=bpy.context.scene.frame_end):
		print('Current frame is '+str(cnt)+' of '+str(end_frame))
		bpy.context.scene.frame_set(cnt)
		bpy.data.scenes["Scene"].frame_current=cnt
		if cnt == (end_frame):
			goNext = check_particles(dimZ, num, goNext)	

		cnt=cnt+1

finishSimulation = time.time()
print("Time to perform simulation: %.3f" % (finishSimulation - finishAddBox))

for obj in name_list:
	delete(obj)
objct = bpy.data.objects
for o in range(np.size(objct)-5):
	obj="p"+str(o)
	[A, B, C]=bpy.data.objects[obj].matrix_world.to_translation()
	out.write("P) "+str(A)+" "+str(B)+" "+str(C)+"\n")
for o in range(np.size(objct)-5):
	obj="p"+str(o)
	[A, B, C]=bpy.data.objects[obj].matrix_world.to_euler()
	out.write("R) "+str(A*360/(2*math.pi))+" "+str(B*360/(2*math.pi))+" "+str(C*360/(2*math.pi))+"\n")
out.write("EOF")
out.close()

if rescale:
	for i in range(num):
		obj="p"+str(i)
		scale(obj,perc=0.95)
		print(i)
export_stl()

finishExport = time.time()
print("Time to export stl file: %.3f" % (finishExport - finishSimulation))
print("Total time: %.3f\n" % (finishExport - startTime))
