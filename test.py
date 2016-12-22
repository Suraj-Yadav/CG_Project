import bpy
import bmesh

bpy.ops.object.delete(use_global=False)

mesh = bmesh.new()

inputFile = open("G:\work\CG_Project\output.txt")

verticesCount = int(inputFile.readline())
vertices = []
edges = []

for i in range(verticesCount):
	coord = tuple(map(float,inputFile.readline().split()))
	vertices.append(mesh.verts.new(coord))

edgeCount = int(inputFile.readline())	
for i in range(edgeCount):
	(u,v) = tuple(map(int,inputFile.readline().split()))
	if (vertices[u],vertices[v]) not in mesh.edges:
		if u>=len(vertices) or v>=len(vertices):
			print(u,v)
		mesh.edges.new((vertices[u],vertices[v]))
#		edges.append((vertices[u],vertices[v]))

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