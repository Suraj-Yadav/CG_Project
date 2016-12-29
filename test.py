import bpy
import bmesh
from typing import List
import sys

inputFile = open(sys.argv[sys.argv.index("--") + 1])

# inputFile = open("G:/work/CG_Project/build/output.txt")

scene = bpy.context.scene
print("scene.objects.active", scene.objects.active)
if scene.objects.active != None:
	bpy.ops.object.delete(use_global=False)

mesh = bmesh.new()	# type: bmesh.BMesh

verticesCount = int(inputFile.readline())	# type: int
vertices = [] # type: List[bmesh.BMVert]

for i in range(verticesCount):
	coord = tuple(map(float,inputFile.readline().split()))
	vertices.append(mesh.verts.new(coord))

edgeCount = int(inputFile.readline())	# type: int

for i in range(edgeCount):
	(u,v) = tuple(map(int,inputFile.readline().split()))
	if (vertices[u],vertices[v]) not in mesh.edges:
		mesh.edges.new((vertices[u],vertices[v]))

faceCount = int(inputFile.readline())	# type: int

for i in range(faceCount):
	(a,b,c) = tuple(map(int,inputFile.readline().split()))
	if (vertices[a], vertices[b], vertices[c]) not in mesh.faces:
		mesh.faces.new((vertices[a], vertices[b], vertices[c]))

me = bpy.data.meshes.new("Mesh")

mesh.to_mesh(me)
mesh.free()

scene = bpy.context.scene
obj = bpy.data.objects.new("Object", me)
scene.objects.link(obj)

bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY')
bpy.ops.object.location_clear(clear_delta = False)


# Select and make active
scene.objects.active = obj
obj.select = True

print("ALL is Well")