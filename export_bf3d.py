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
import numpy as np
from bpy_extras.io_utils import axis_conversion
from bpy.props import *
from mathutils import Vector, Quaternion, Matrix
from . import struct_bf3d

#TODO 

# export the hierarchy to the model file if no skeleton exists

# animation export

HEAD = 8 #4(int = chunktype) + 4 (int = chunksize)

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
	WriteFloat(file, vec[0])
	WriteFloat(file, vec[1])
	WriteFloat(file, vec[2])

def WriteQuaternion(file, quat):
	WriteFloat(file, quat[0])
	WriteFloat(file, quat[1])
	WriteFloat(file, quat[2])
	WriteFloat(file, quat[3])
	
def WriteMatrix(file, mat):
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
	
def WriteAnimationHeader(file, size, header):
	WriteInt(file, 513) #chunktype
	WriteInt(file, size) #chunksize

	WriteString(file, header.name)
	WriteString(file, header.hieraName)
	WriteInt(file, header.numFrames)
	WriteInt(file, header.frameRate)

def WriteTimeCodedAnimationChannel(file, channel):
	WriteInt(file, 514) #chunktype
	size = 6 + (len(channel.timeCodedKeys) * channel.vectorLen) * 4
	WriteInt(file, size) #chunksize
	
	WriteInt(file, channel.vectorLen)
	WriteInt(file, channel.type)
	WriteInt(file, channel.pivot)

	if channel.vectorLen == 1:
		for f in channel.timeCodedKeys:
			WriteFloat(file, f)
	elif channel.vectorLen == 4:
		for quat in channel.timeCodedKeys:
			WriteQuaternion(file, quat)

def WriteAnimation(file, animation):
	print("\n### NEW ANIMATION: ###")
	WriteInt(file, 512) #chunktype
	
	headerSize = len(header.name) + len(header.hieraName) + 8
	channelsSize = 0
	for channel in animation.channels:
		channelsSize += HEAD + 6 + (len(channel.timeCodedKeys) * channel.vectorLen) * 4
	size = HEAD + headerSize + channelsSize		
	
	WriteInt(file, size) #chunksize
	
	WriteAnimationHeader(file, headerSize, animation.header)
	print("Header")
	for channel in animation.channels:
		WriteAnimationChannel(file, channel)
		print("Channel")
			
#######################################################################################
# Box
#######################################################################################

def getBoxChunkSize(box):
	return 24

def WriteBox(file, box):
	print("\n### NEW BOX: ###")
	WriteInt(file, 1024) #chunktype
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
	WriteInt(file, 3) #chunktype
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
	WriteInt(file, 4) #chunktype
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
	WriteInt(file, 5) #chunktype
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
	WriteInt(file, 6) #chunktype
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
	WriteInt(file, 7) #chunktype
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
	WriteInt(file, 2) #chunktype
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
	WriteInt(file, 1) #chunktype

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
	if not model.bVolume == None:
		size += getBoxChunkSize(model.bVolume)
	for mesh in model.meshes:
		size += getMeshChunkSize(mesh)
	return size

def WriteModel(file, model):
	print("\n### NEW MODEL: ###")
	WriteInt(file, 0) #chunktype
	WriteInt(file, getModelChunkSize(model)) #chunksize

	print(model.hieraName)
	WriteString(file, model.hieraName)
	if not model.bVolume == None:
		WriteBox(file, model.bVolume)
	for mesh in model.meshes:
		WriteMesh(file, mesh)
		
#######################################################################################
# Main Export
#######################################################################################

def MainExport(givenfilepath, self, context, EXPORT_MODE = 'M'):
	#axis_conversion(from_forward='Y', from_up='Z', to_forward='Z', to_up='Y')
	#print("Run Export")
	Hierarchy = struct_bf3d.Hierarchy()
	Hierarchy.pivots = []
	amtName = ""
	modelName = ""
 
	roottransform = struct_bf3d.HierarchyPivot()
	roottransform.name = "ROOTTRANSFORM"
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
 
	if EXPORT_MODE == 'M' or EXPORT_MODE == 'H':
		modelName = (os.path.splitext(os.path.basename(givenfilepath))[0]).upper()
		print(modelName)
 
		Model = struct_bf3d.Model()
		Model.name = modelName
		Model.hieraName = amtName
		Model.meshes = []
 
		for mesh_ob in objList: 
			if mesh_ob.name == "BOUNDINGBOX":
				Box = struct_bf3d.Box()
				Box.center = mesh_ob.location
				box_mesh = mesh_ob.to_mesh(bpy.context.scene, False, 'PREVIEW', calc_tessface = True)
				Box.extend = Vector((box_mesh.vertices[0].co.x * 2, box_mesh.vertices[0].co.y * 2, box_mesh.vertices[0].co.z))
 
				if not EXPORT_MODE == 'H':
					Model.bVolume = Box
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
					pivot.isBone = 0
					pivot.matrix = mesh_ob.matrix_world
					if not mesh_ob.parent_bone == "":
						ids = [index for index, pivot in enumerate(Hierarchy.pivots) if pivot.name == mesh_ob.parent_bone] #return an array of indices (in this case only one value)
						pivot.parent = ids[0]
					pivot.isBone = 0
					Mesh.header.parentPivot = len(Hierarchy.pivots)
					Hierarchy.pivots.append(pivot)

				Model.meshes.append(Mesh)

		if EXPORT_MODE == 'M':
			sknFile = open(givenfilepath, "wb")
			WriteModel(sknFile, Model)
			sknFile.close()

	Hierarchy.header.pivotCount = len(Hierarchy.pivots)

	if EXPORT_MODE == 'H':
		sklFile = open(givenfilepath.replace(modelName.lower(), amtName.lower()), "wb")
		Hierarchy.header.name = amtName
		WriteHierarchy(sklFile, Hierarchy) 
		sklFile.close()