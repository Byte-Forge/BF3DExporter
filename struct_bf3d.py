#Written by Michael Schnabel
#Last Modification 08.02.2016
#Structs of the BF3D Format 
from mathutils import Vector, Quaternion

class Struct:
    def __init__ (self, *argv, **argd):
        if len(argd):
            # Update by dictionary
            self.__dict__.update (argd)
        else:
            # Update by position
            attrs = filter (lambda x: x[0:2] != "__", dir(self))
            for n in range(len(argv)):
                setattr(self, attrs[n], argv[n])
			
#######################################################################################
# Basic Structs
#######################################################################################
			
class RGBA(Struct):
    r = 0
    g = 0
    b = 0
    a = 0
	
#######################################################################################
# Hierarchy
#######################################################################################

#chunk 257
class HierarchyHeader(Struct):
    name = ""
    pivotCount = 0
    centerPos = Vector((0.0, 0.0 ,0.0))

#chunk 258
class HierarchyPivot(Struct):
    name = ""
    parentID = -1
    isBone = 1 #default 1
    position = Vector((0.0, 0.0 ,0.0))
    rotation = Quaternion((1.0, 0.0, 0.0, 0.0))

# chunk 256
class Hierarchy(Struct):
    header = HierarchyHeader()
    pivots = []
	
#######################################################################################
# Animation
#######################################################################################

#chunk 513
class AnimationHeader(Struct):
    name = ""
    hieraName = ""
    numFrames = 0
    frameRate = 0
	
class TimeCodedAnimationKey(Struct):
    frame = 0
    value = 0
	
#chunk 514
class TimeCodedAnimationChannel(Struct):
    vectorLen = 0
    type = 0
    pivot = 0 
    timeCodedKeys = []
	
#chunk 512
class Animation(Struct):
    header = AnimationHeader()
    channels = [] 
	
#######################################################################################
# Model
#######################################################################################

#chunk 0
class Model(Struct):
    hieraName = "" # is empty
    meshes = []
    bVolume = None
	
#######################################################################################
# Box
#######################################################################################	

#chunk 1024
class Box(Struct): 
    center = Vector((0.0, 0.0 ,0.0))
    extend = Vector((0.0, 0.0 ,0.0))
	
#######################################################################################
# VertexInfluences
#######################################################################################

#chunk 7
class MeshVertexInfluences(Struct):
    boneIdx = 0
    boneInf = 0.0
	
#######################################################################################
# Mesh
#######################################################################################	

#chunk 2
class MeshHeader(Struct):
    type = 0
    # 0   -> normal mesh
	# 1   -> normal mesh - two sided
    # 2   -> normal mesh - camera oriented
    # 128 -> skin
	# 129 -> skin - two sided
   
    meshName = ""
    materialID = 0
    parentPivot = 0
    faceCount = 0
    vertCount = 0

#chunk 1
class Mesh(Struct):
    header = MeshHeader()
    verts = []
    normals = []
    faces = []
    uvCoords = []
    vertInfs = []