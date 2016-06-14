#Written by Michael Schnabel
#Last Modification 08.02.2016
#Exports the BF3D Format 
import bpy
import operator
import struct
import os
import math
import sys
import bmesh
from bpy_extras.io_utils import axis_conversion
from bpy.props import *
from mathutils import Vector, Quaternion, Matrix
from . import struct_bf3d

#TODO 

#compute a minimal bounding sphere

version = 1.0

HEAD = 8 #4(int = chunktype) + 4 (int = chunksize)

#matrix for axis conversion from z up to y up
#for animations we have to do this manually !!!
global_matrix = axis_conversion(from_forward='Y', from_up='Z', to_forward='-Z', to_up='Y').to_4x4()

#######################################################################################
# Basic Methods
#######################################################################################

def getStringSize(string):
	return len(string) + 1 #binary 0

def WriteString(file, string):
	file.write(bytes(string, 'UTF-8'))
	#write binary 0 to file
	file.write(struct.pack('B', 0))
		
def WriteRGBA(file, rgba):
	file.write(struct.pack("B", int(rgba.r)))
	file.write(struct.pack("B", int(rgba.g)))
	file.write(struct.pack("B", int(rgba.b)))
	file.write(struct.pack("B", int(rgba.a)))

def WriteInt(file, num):
	file.write(struct.pack("<i", num))
	
def WriteIntArray(file, array):
	for a in array:
		WriteLong(file, a)

def WriteFloat(file, num):
	file.write(struct.pack("<f", num))

def WriteUnsignedByte(file, num):
	file.write(struct.pack("<B", num))

def WriteVector(file, vec):
	vec = global_matrix * vec
	WriteFloat(file, vec[0])
	WriteFloat(file, vec[1])
	WriteFloat(file, vec[2])

def WriteQuaternion(file, quat):
	WriteFloat(file, quat[0])
	WriteFloat(file, quat[1])
	WriteFloat(file, quat[2])
	WriteFloat(file, quat[3])
	
def WriteMatrix(file, mat):
	mat = global_matrix * mat * global_matrix.inverted()
	WriteQuaternion(file, mat[0])
	WriteQuaternion(file, mat[1])
	WriteQuaternion(file, mat[2])
	WriteQuaternion(file, mat[3])
	
#######################################################################################
# Triangulate
#######################################################################################	

def triangulate(mesh):
	import bmesh
	bm = bmesh.new()
	bm.from_mesh(mesh)
	bmesh.ops.triangulate(bm, faces = bm.faces)
	bm.to_mesh(mesh)
	bm.free()
	
#######################################################################################
# Mesh Sphere
#######################################################################################

def calcModelSphere(model):
	objList = [object for object in bpy.context.scene.objects if object.type == 'MESH']
	verts = []
	for mesh_ob in objList: 
		if (mesh_ob.name == 'BOUNDINGBOX'):
			continue
		for b in mesh_ob.bound_box:
			verts.append(mesh_ob.matrix_world * Vector(b))
	
	# get the point with the biggest distance to x and store it in y
	x = verts[0]
	y = verts[1]
	z = verts[2]
	
	dist = ((y.x- x.x)**2 + (y.y - x.y)**2 + (y.z - x.z)**2)**(1/2)
	for v in verts:
		curr_dist = ((v.x - x.x)**2 + (v.y - x.y)**2 + (v.z - x.z)**2)**(1/2)
		if (curr_dist > dist):
			dist = curr_dist
			y = v
					
	#get the point with the biggest distance to y and store it in z
	dist = ((z.x - y.x)**2 + (z.y - y.y)**2 + (z.z - y.z)**2)**(1/2)
	for v in verts:
		curr_dist = ((v.x - y.x)**2 + (v.y - y.y)**2 + (v.z - y.z)**2)**(1/2)
		if (curr_dist > dist):
			dist = curr_dist
			z = v
			 
	# the center of the sphere is between y and z
	y_z = ((z - y)/2)
	m = (y + z)/2
	radius = y_z.length

	
	#test if any of the vertices is outside the sphere (if so update the sphere)
	for v in verts:
		curr_dist = ((v.x - m.x)**2 + (v.y - m.y)**2 + (v.z - m.z)**2)**(1/2)
		if curr_dist > radius:
			delta = (curr_dist - radius)/2
			radius += delta
			m += (Vector(v - m)).normalized() * delta
	s = struct_bf3d.Sphere()
	s.center = m
	s.extend = radius
	model.bSphere = s
	
	#for testing
	createSphere(radius, m.x, m.y, m.z)
	
#######################################################################################
# create Sphere (just for testing)
#######################################################################################

def createSphere(radius, x, y, z):
	bpy.ops.mesh.primitive_uv_sphere_add(size=radius, location=(x, y, z))
	objList = [object for object in bpy.context.scene.objects if object.type == 'MESH']
	for mesh_ob in objList: 
		mesh_ob.draw_type = 'WIRE'
	
#######################################################################################
# BF3D
#######################################################################################

def WriteBF3D(file, name):
	file.write(bytes("BF3D", 'UTF-8'))
	WriteInt(file, 0) #chunktype
	WriteInt(file, 4) #chunksize
	WriteFloat(file, version)
	
#######################################################################################
# Hierarchy
#######################################################################################

def getHierarchyHeaderChunkSize(header):
	return 16 + getStringSize(header.name)

def WriteHierarchyHeader(file, header):
	WriteInt(file, 257) #chunktype
	WriteInt(file, getHierarchyHeaderChunkSize(header)) #chunksize
	
	WriteString(file, header.name)
	WriteInt(file, header.pivotCount)
	WriteVector(file, header.centerPos)
	
def getPivotsChunkSize(pivots):
	size = 0
	for pivot in pivots:
		size += 69 + getStringSize(pivot.name)
	return size

def WritePivots(file, pivots):
	WriteInt(file, 258) #chunktype
	WriteInt(file, getPivotsChunkSize(pivots)) #chunksize
	
	for pivot in pivots:
		WriteString(file, pivot.name)
		WriteInt(file, pivot.parent)
		WriteUnsignedByte(file, pivot.isBone)
		WriteMatrix(file, pivot.matrix)

def WriteHierarchy(file, hierarchy):
	print("\n### NEW HIERARCHY: ###")
	WriteInt(file, 256) #chunktype
	
	size = HEAD + getHierarchyHeaderChunkSize(hierarchy.header) + HEAD + getPivotsChunkSize(hierarchy.pivots) 

	WriteInt(file, size) #chunksize
	
	WriteHierarchyHeader(file, hierarchy.header)
	print("Header")
	WritePivots(file, hierarchy.pivots)
	print("Pivots")

#######################################################################################
# Animation
#######################################################################################

def getAnimationHeaderChunkSize(header):
	return 8 + getStringSize(header.name) + getStringSize(header.hieraName)
	
def WriteAnimationHeader(file, header):
	WriteInt(file, 513) #chunktype
	WriteInt(file, getAnimationHeaderChunkSize(header)) #chunksize

	WriteString(file, header.name)
	WriteString(file, header.hieraName)
	WriteFloat(file, header.frameRate)
	WriteInt(file, header.numFrames)
	
def getTimeCodedAnimationChannelSize(channel):
	size = 12
	size += len(channel.timeCodedKeys) * 8
	return size

def WriteTimeCodedAnimationChannel(file, channel):
	WriteInt(file, 514) #chunktype
	WriteInt(file, getTimeCodedAnimationChannelSize(channel)) #chunksize

	WriteInt(file, channel.pivot)
	WriteInt(file, channel.extrapolation)
	WriteInt(file, channel.type)
	
	for key in channel.timeCodedKeys:
		WriteInt(file, int(key.frame))
		WriteFloat(file, key.value)

def WriteAnimation(file, animation):
	print("\n### NEW ANIMATION: ###")
	WriteInt(file, 512) #chunktype
	channelsSize = 0
	for channel in animation.channels:
		channelsSize += HEAD + getTimeCodedAnimationChannelSize(channel)
	WriteInt(file, HEAD + getAnimationHeaderChunkSize(animation.header) + channelsSize) #chunksize
	
	WriteAnimationHeader(file, animation.header)
	print("Header")
	for channel in animation.channels:
		WriteTimeCodedAnimationChannel(file, channel)
		
#######################################################################################
# Sphere
#######################################################################################

def getSphereChunkSize(sphere):
	return 16

def WriteSphere(file, sphere):
	print("\n### NEW SPHERE: ###")
	WriteInt(file, 193) #chunktype
	WriteInt(file, getSphereChunkSize(sphere)) #chunksize
	
	WriteVector(file, sphere.center)
	WriteFloat(file, sphere.radius)
			
#######################################################################################
# Box
#######################################################################################

def getBoxChunkSize(box):
	return 24

def WriteBox(file, box):
	print("\n### NEW BOX: ###")
	WriteInt(file, 192) #chunktype
	WriteInt(file, getBoxChunkSize(box)) #chunksize
	
	WriteVector(file, box.center)
	WriteVector(file, box.extend)
	
#######################################################################################
# Vertices
#######################################################################################

def getMeshVerticesChunkSize(vertices):
	size = len(vertices) * 12
	return size

def WriteMeshVerticesArray(file, vertices):
	WriteInt(file, 131) #chunktype
	WriteInt(file, getMeshVerticesChunkSize(vertices)) #chunksize
	
	for vert in vertices:
		WriteVector(file, vert)

#######################################################################################
# Normals
#######################################################################################

def getMeshNormalsArrayChunkSize(normals):
	size = len(normals) * 12
	return size
	
def WriteMeshNormalsArray(file, normals):
	WriteInt(file, 132) #chunktype
	WriteInt(file, getMeshNormalsArrayChunkSize(normals)) #chunksize
	
	for norm in normals:
		WriteVector(file, norm)
	
#######################################################################################
# Faces
#######################################################################################	

def getMeshFaceArrayChunkSize(faces):
	size = len(faces) * 12
	return size

def WriteMeshFaceArray(file, faces):
	WriteInt(file, 133) #chunktype
	WriteInt(file, getMeshFaceArrayChunkSize(faces)) #chunksize
	
	for face in faces:
		WriteInt(file, face[0])
		WriteInt(file, face[1])
		WriteInt(file, face[2])
		
#######################################################################################
# uvCoords
#######################################################################################	

def getMeshUVCoordsChunkSize(uvCoords):
	return len(uvCoords) * 8

def WriteMeshUVCoords(file, uvCoords):
	WriteInt(file, 134) #chunktype
	WriteInt(file, getMeshUVCoordsChunkSize(uvCoords)) #chunksize
	
	for uv in uvCoords:
		WriteFloat(file, uv[0])
		WriteFloat(file, uv[1])
		
#######################################################################################
# VertexInfluences
#######################################################################################	
		
def getMeshVertexInfluencesChunkSize(influences):
	size = len(influences) * 8
	return size

def WriteMeshVertexInfluences(file, influences):
	WriteInt(file, 135) #chunktype
	WriteInt(file, getMeshVertexInfluencesChunkSize(influences)) #chunksize

	for inf in influences:
		WriteInt(file, inf.boneIdx)
		WriteInt(file, int(inf.boneInf * 100))
		
#######################################################################################
# Mesh
#######################################################################################	

def getMeshHeaderChunkSize(header):
	size = 17 + getStringSize(header.meshName)
	return size

def WriteMeshHeader(file, header): 
	WriteInt(file, 130) #chunktype
	WriteInt(file, getMeshHeaderChunkSize(header)) #chunksize

	WriteUnsignedByte(file, header.type)
	WriteString(file, header.meshName)
	WriteInt(file, header.materialID)
	WriteInt(file, header.parentPivot)
	WriteInt(file, header.faceCount)
	WriteInt(file, header.vertCount)
	
def getMeshChunkSize(mesh):
	size = HEAD + getMeshHeaderChunkSize(mesh.header)
	size += HEAD + getMeshVerticesChunkSize(mesh.verts)
	size += HEAD + getMeshNormalsArrayChunkSize(mesh.normals)
	size += HEAD + getMeshFaceArrayChunkSize(mesh.faces)
	size += HEAD + getMeshUVCoordsChunkSize(mesh.uvCoords)
	if len(mesh.vertInfs) > 0:
		size += HEAD + getMeshVertexInfluencesChunkSize(mesh.vertInfs)
	return size
	
def WriteMesh(file, mesh):
	print("\n### NEW MESH: ###")
	WriteInt(file, 129) #chunktype

	WriteInt(file, getMeshChunkSize(mesh)) #chunksize
	
	WriteMeshHeader(file, mesh.header)
	print(mesh.header.meshName)
	#print("Header")
	WriteMeshVerticesArray(file, mesh.verts)
	#print("Vertices")
	WriteMeshNormalsArray(file, mesh.normals)
	#print("Normals")
	WriteMeshFaceArray(file, mesh.faces)
	#print("Faces")
	WriteMeshUVCoords(file, mesh.uvCoords)
	#print("uvCoords")
	if len(mesh.vertInfs) > 0:
		WriteMeshVertexInfluences(file, mesh.vertInfs) 
		#print("Vertex Influences")
		
#######################################################################################
# Model
#######################################################################################

def getModelChunkSize(model):
	size = getStringSize(model.hieraName)
	if not model.bBox == None:
		size += getBoxChunkSize(model.bBox)
	if not model.bSphere == None:
		size += getSphereChunkSize(model.bSphere)
	for mesh in model.meshes:
		size += getMeshChunkSize(mesh)
	return size

def WriteModel(file, model):
	print("\n### NEW MODEL: ###")
	WriteInt(file, 128) #chunktype
	WriteInt(file, getModelChunkSize(model)) #chunksize

	calcModelSphere(model)
	print(model.hieraName)
	WriteString(file, model.hieraName)
	if not model.bBox == None:
		WriteBox(file, model.bBox)
	if not model.bSphere == None:
		WriteSphere(file, model.bSphere)
	for mesh in model.meshes:
		WriteMesh(file, mesh)
		
#######################################################################################
# Main Export
#######################################################################################

def MainExport(givenfilepath, self, context, EXPORT_MODE = 'M'):
	#print("Run Export")
	fileName = os.path.splitext(os.path.basename(givenfilepath))[0]
	Hierarchy = struct_bf3d.Hierarchy()
	Animation = struct_bf3d.Animation()
	Hierarchy.pivots = []
	amtName = ""
	modelName = ""
	pivotsList = []
 
	roottransform = struct_bf3d.HierarchyPivot()
	roottransform.name = "ROOTTRANSFORM"
	pivotsList.append(roottransform.name)
	roottransform.matrix = Matrix()
	roottransform.parent = -1
	Hierarchy.pivots.append(roottransform)
	
	#switch to object mode
	if bpy.ops.object.mode_set.poll():
		bpy.ops.object.mode_set(mode='OBJECT')

	# Get all the armatures in the scene.
	rigList = [object for object in bpy.context.scene.objects if object.type == 'ARMATURE']

	if len(rigList) == 1:
		rig = rigList[0]
		amtName = rig.name
		
		for bone in rig.pose.bones:
			pivot = struct_bf3d.HierarchyPivot()
			pivot.name = bone.name
			pivotsList.append(bone.name)
			pivot.parent = 0
			pivot.matrix = bone.matrix_basis.copy()
			if not bone.parent == None:
				ids = [index for index, pivot in enumerate(Hierarchy.pivots) if pivot.name == bone.parent.name] #return an array of indices (in this case only one value)
				pivot.parent = ids[0]
			Hierarchy.pivots.append(pivot)
			
	if len(rigList) > 1:
		context.report({'ERROR'}, "only one armature allowed!")
		print("Error: only one armature allowed!") 

	objList = []
	# Get all the mesh objects in the scene.
	objList = [object for object in bpy.context.scene.objects if object.type == 'MESH']
 
	modelName = fileName

	Model = struct_bf3d.Model()
	Model.name = modelName
	Model.hieraName = amtName
	Model.meshes = []

	for mesh_ob in objList: 
		if mesh_ob.name == "BOUNDINGBOX":
			Box = struct_bf3d.Box()
			Box.center = (mesh_ob.matrix_world * Vector(mesh_ob.bound_box[0]) + mesh_ob.matrix_world * Vector(mesh_ob.bound_box[6])) / 2.0
			Box.extend = Box.center - mesh_ob.matrix_world * Vector(mesh_ob.bound_box[0])
			Model.bBox = Box
		else:
			Mesh = struct_bf3d.Mesh()
			Mesh.header = struct_bf3d.MeshHeader()
			Mesh.verts = []
			Mesh.normals = [] 
			Mesh.faces = []
			Mesh.uvCoords = []
			Mesh.vertInfs = []

			Mesh.header.meshName = mesh_ob.name
			mesh = mesh_ob.to_mesh(bpy.context.scene, False, 'PREVIEW', calc_tessface = True)

			triangulate(mesh)

			Mesh.header.vertCount = len(mesh.vertices)
  
			group_lookup = {g.index: g.name for g in mesh_ob.vertex_groups}
			groups = {name: [] for name in group_lookup.values()}

			for face in mesh.polygons:
				Mesh.faces.append((face.vertices[0], face.vertices[1], face.vertices[2]))

			Mesh.header.faceCount = len(Mesh.faces)

			for v in mesh.vertices:
				Mesh.verts.append(v.co.xyz)
				Mesh.normals.append(v.normal)
				Mesh.uvCoords.append((0.0, 0.0)) #just to fill the array 

				#vertex influences
				vertInf = struct_bf3d.MeshVertexInfluences()
				if len(v.groups) == 1:
					#has to be this complicated, otherwise the vertex groups would be corrupted
					ids = [index for index, pivot in enumerate(Hierarchy.pivots) if pivot.name == mesh_ob.vertex_groups[v.groups[0].group].name] #return an array of indices (in this case only one value)
					if len(ids) > 0:
						vertInf.boneIdx = ids[0]
					vertInf.boneInf = v.groups[0].weight
					Mesh.vertInfs.append(vertInf)
				elif len(v.groups) == 2:
					#has to be this complicated, otherwise the vertex groups would be corrupted
					ids = [index for index, pivot in enumerate(Hierarchy.pivots) if pivot.name == mesh_ob.vertex_groups[v.groups[0].group].name] #return an array of indices (in this case only one value)
					if len(ids) > 0:
						vertInf.boneIdx = ids[0]
					vertInf.boneInf = v.groups[0].weight
					#has to be this complicated, otherwise the vertex groups would be corrupted
					ids = [index for index, pivot in enumerate(Hierarchy.pivots) if pivot.name == mesh_ob.vertex_groups[v.groups[1].group].name] #return an array of indices (in this case only one value)
					if len(ids) > 0:
						vertInf.xtraIdx = ids[0]
					vertInf.xtraInf = v.groups[1].weight
					Mesh.vertInfs.append(vertInf)
				elif len(v.groups) > 2: 
					context.report({'ERROR'}, "max 2 bone influences per vertex supported!")
					print("Error: max 2 bone influences per vertex supported!")

			#uv coords
			bm = bmesh.new()
			bm.from_mesh(mesh)

			uv_layer = bm.loops.layers.uv.verify() 

			index = 0
			for f in bm.faces:
				#test if we need this 1- at all meshes
				Mesh.uvCoords[Mesh.faces[index][0]] = (f.loops[0][uv_layer].uv[0], 1 - f.loops[0][uv_layer].uv[1])
				Mesh.uvCoords[Mesh.faces[index][1]] = (f.loops[1][uv_layer].uv[0], 1 - f.loops[1][uv_layer].uv[1])
				Mesh.uvCoords[Mesh.faces[index][2]] = (f.loops[2][uv_layer].uv[0], 1 - f.loops[2][uv_layer].uv[1])
				index+=1   
			del bm
			
			if len(mesh_ob.vertex_groups) > 0:
				Mesh.header.type = 128 #type skin
			else:
				Mesh.header.type = 0 #type normal mesh
				pivot = struct_bf3d.HierarchyPivot()
				pivot.name = mesh_ob.name
				pivotsList.append(mesh_ob.name)
				pivot.isBone = 0
				pivot.matrix = mesh_ob.matrix_basis
				if not mesh_ob.parent_bone == "":
					ids = [index for index, pivot in enumerate(Hierarchy.pivots) if pivot.name == mesh_ob.parent_bone] #return an array of indices (in this case only one value)
					pivot.parent = ids[0]
				elif not mesh_ob.parent == None:
					ids = [index for index, pivot in enumerate(Hierarchy.pivots) if pivot.name == mesh_ob.parent.name] #return an array of indices (in this case only one value)
					pivot.parent = ids[0]
				pivot.isBone = 0
				Mesh.header.parentPivot = len(Hierarchy.pivots)
				Hierarchy.pivots.append(pivot)

			Model.meshes.append(Mesh)

	if EXPORT_MODE == 'M':
		sknFile = open(givenfilepath, "wb")
		WriteBF3D(sknFile, fileName)
		WriteModel(sknFile, Model)
		sknFile.close()

	Hierarchy.header.pivotCount = len(Hierarchy.pivots)
	
	if EXPORT_MODE == 'A':
		if len(rigList) == 1: #could also be 0?
			Animation.header = struct_bf3d.AnimationHeader()
			Animation.header.hieraName = amtName
			Animation.header.frameRate = bpy.data.scenes["Scene"].render.fps
			Animation.header.numFrames = bpy.data.scenes["Scene"].frame_end - bpy.data.scenes["Scene"].frame_start
			Animation.channels = []
			for obj in bpy.data.objects:
				if obj.animation_data == None:
					continue
				action = obj.animation_data.action
				frame_begin, frame_end = [int(x) for x in action.frame_range]

				for fcu in action.fcurves:
					channel = struct_bf3d.TimeCodedAnimationChannel()
					if(fcu.extrapolation == "CONSTANT"):
						channel.extrapolation = 1
					elif(fcu.extrapolation == "BEIZIER"):
						channel.extrapolation = 2
					channel.type = fcu.array_index 

					if (fcu.data_path.endswith("location")):
						channel.type += 0
					elif (fcu.data_path.endswith("rotation_quaternion")):
						channel.type += 3
					else:
						print("ERROR!: that type of data_path is not supported yet!")
						print(fcu.data_path)
						continue
					channel.timeCodedKeys = []
					try:
						pivotName = fcu.data_path.split('"')[1]
					except:
						pivotName = obj.name
					channel.pivot = pivotsList.index(pivotName)
					
					#axis conversion is applied here
					if channel.type == 1:
						channel.type = 2
					elif channel.type == 2:
						channel.type = 1
					elif channel.type == 5:
						channel.type = 6
					elif channel.type == 6:
						channel.type = 5

					for keyframe in fcu.keyframe_points:
						key = struct_bf3d.TimeCodedAnimationKey()
						key.frame = keyframe.co.x

						if channel.type == 0:
							key.value = keyframe.co.y - Hierarchy.pivots[channel.pivot].matrix[0][3]
						elif channel.type == 1:
							key.value = keyframe.co.y - Hierarchy.pivots[channel.pivot].matrix[2][3]
						elif channel.type == 2:
							key.value = -(keyframe.co.y - Hierarchy.pivots[channel.pivot].matrix[1][3])
						
						elif channel.type == 3:
							key.value = keyframe.co.y
						elif channel.type == 4:
							key.value = -keyframe.co.y
						elif channel.type == 5:
							key.value = -keyframe.co.y
						elif channel.type == 6:
							key.value = keyframe.co.y
						else:
							print("invalid animation channel type")
						channel.timeCodedKeys.append(key)
					Animation.channels.append(channel)

	if EXPORT_MODE == 'H':
		sklFile = open(givenfilepath.replace(fileName, amtName), "wb")
		Hierarchy.header.name = amtName
		WriteBF3D(sklFile, amtName)
		WriteHierarchy(sklFile, Hierarchy) 
		sklFile.close()
		
	if EXPORT_MODE == 'A':
		aniFile = open(givenfilepath, "wb")
		WriteBF3D(aniFile, fileName)
		WriteAnimation(aniFile, Animation)
		aniFile.close()