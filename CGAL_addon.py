bl_info = {'name': 'CGAL Integration', 'category': 'All'}

import bpy
import subprocess
from bpy.types import Operator
from bpy.types import Panel
from bpy.types import Context
import bmesh
from typing import List
from os import path

addon_keymaps = []


def showOnlyThis(objName: str):
	if objName in bpy.data.objects:
		for obj in bpy.data.objects:
			obj.select = False
			obj.hide = True
		bpy.data.objects[objName].hide = False
		bpy.context.scene.objects.active = bpy.data.objects[objName]
		bpy.data.objects[objName].select = True
		return True
	return False


def selectInMesh(obj, selectionType: str, indices: bytes):
	if obj.mode == 'EDIT':
		bpy.ops.object.editmode_toggle()

	index = indices.decode().split()  # type: List[str]

	elements = getattr(obj.data, selectionType)

	for i in range(len(elements)):
		elements[i].select = False

	for iStr in index:  # type: str
		i = -1
		if iStr.isdigit():
			i = int(iStr)
		if i < 0 or i >= len(elements):
			print("'" + iStr + "'", "is an Invalid Index")
		else:
			elements[i].select = True

	bpy.ops.object.editmode_toggle()


def loadModelFromFile(outputFile: str, context):
	scene = context.scene

	for obj in bpy.data.objects:
		obj.hide = False
		obj.select = False

	for obj in bpy.data.objects:
		obj.select = True
		bpy.ops.object.delete(use_global=False)

	edgeMesh = bmesh.new()  # type: bmesh.BMesh
	faceMesh = bmesh.new()  # type: bmesh.BMesh

	vertices = []  # type: List[Tuple[float]]
	edges = []  # type: List[Tuple[int]]
	faces = []  # type: List[Tuple[int]]

	inputFile = open(outputFile)

	verticesCount = int(inputFile.readline())  # type: int
	for i in range(verticesCount):
		vertices.append(list(map(float, inputFile.readline().split())))

	edgesCount = int(inputFile.readline())  # type: int
	for i in range(edgesCount):
		edges.append(list(map(int, inputFile.readline().split())))

	facesCount = int(inputFile.readline())  # type: int
	for i in range(facesCount):
		faces.append(list(map(int, inputFile.readline().split())))

	if verticesCount > 0:
		me = bpy.data.meshes.new("Vertex_Mesh")
		vertexMesh = bmesh.new()  # type: bmesh.BMesh
		vertexHandle = [vertexMesh.verts.new(i) for i in vertices]
		vertexMesh.to_mesh(me)
		vertexMesh.free()
		obj = bpy.data.objects.new("Model_Vertex", me)
		scene.objects.link(obj)
		obj.data.show_edge_seams = True
		print("Vertex Model Generated")

	if edgesCount > 0:
		me = bpy.data.meshes.new("Edges_Mesh")
		edgeMesh = bmesh.new()  # type: bmesh.BMesh
		vertexHandle = [edgeMesh.verts.new(i) for i in vertices]
		for edge in edges:
			edgeMesh.edges.new((vertexHandle[edge[0]], vertexHandle[edge[1]]))

		edgeMesh.to_mesh(me)
		edgeMesh.free()
		obj = bpy.data.objects.new("Model_Edge", me)
		scene.objects.link(obj)
		obj.data.show_edge_seams = True
		# obj.show_x_ray = True
		print("Edge Model Generated")

	if facesCount > 0:
		me = bpy.data.meshes.new("Face_Mesh")
		faceMesh = bmesh.new()  # type: bmesh.BMesh
		vertexHandle = [faceMesh.verts.new(i) for i in vertices]
		for face in faces:
			faceMesh.faces.new((vertexHandle[face[0]], vertexHandle[face[1]], vertexHandle[face[2]]))
		# for face in faces:
		# 	a = faceMesh.verts.new(vertices[face[0]])
		# 	b = faceMesh.verts.new(vertices[face[1]])
		# 	c = faceMesh.verts.new(vertices[face[2]])
		# 	faceMesh.faces.new((a, b, c))

		for edge in edges:
			meshEdge = faceMesh.edges.get((vertexHandle[edge[0]], vertexHandle[edge[1]]))
			if meshEdge is None:
				meshEdge = faceMesh.edges.new((vertexHandle[edge[0]], vertexHandle[edge[1]]))
			meshEdge.seam = True

		faceMesh.to_mesh(me)
		faceMesh.free()
		obj = bpy.data.objects.new("Model_Face", me)
		scene.objects.link(obj)
		obj.data.show_edge_seams = True
		print("Face Model Generated")

	return {"FINISHED"}


class vertexSelector(Operator):
	bl_idname = 'cgal.select_given_vertices'
	bl_label = 'Select Vertices'
	bl_options = {"REGISTER", "UNDO"}

	def execute(self, context: Context):
		obj = context.active_object
		if obj == None:
			self.report({'ERROR'}, "No Active Object")
			return {"CANCELLED"}
		selectInMesh(obj, 'vertices', context.scene.vertsToSelect)
		return {"FINISHED"}


class edgeSelector(Operator):
	bl_idname = 'cgal.select_given_edges'
	bl_label = 'Select Edges'
	bl_options = {"REGISTER", "UNDO"}

	def execute(self, context: Context):
		obj = context.active_object
		if obj == None:
			self.report({'ERROR'}, "No Active Object")
			return {"CANCELLED"}
		selectInMesh(obj, 'edges', context.scene.edgesToSelect)
		return {"FINISHED"}


class faceSelector(Operator):
	bl_idname = 'cgal.select_given_faces'
	bl_label = 'Select Faces'
	bl_options = {"REGISTER", "UNDO"}

	def execute(self, context: Context):
		obj = context.active_object
		if obj == None:
			self.report({'ERROR'}, "No Active Object")
			return {"CANCELLED"}
		selectInMesh(obj, 'polygons', context.scene.facesToSelect)
		return {"FINISHED"}


class vertexModel(Operator):
	bl_idname = 'cgal.show_cgal_vertex_model'
	bl_label = 'Show Vertex Model'
	bl_options = {"REGISTER", "UNDO"}

	def execute(self, context: Context):
		if showOnlyThis('Model_Vertex'):
			return {"FINISHED"}
		else:
			self.report({'ERROR'}, "No Vertex Data Found")
			return {"CANCELLED"}


class edgeModel(Operator):
	bl_idname = 'cgal.show_cgal_edge_model'
	bl_label = 'Show Edge Model'
	bl_options = {"REGISTER", "UNDO"}

	def execute(self, context: Context):
		if showOnlyThis('Model_Edge'):
			return {"FINISHED"}
		else:
			self.report({'ERROR'}, "No Edge Data Found")
			return {"CANCELLED"}


class faceModel(Operator):
	bl_idname = 'cgal.show_cgal_face_model'
	bl_label = 'Show Face Model'
	bl_options = {"REGISTER", "UNDO"}

	def execute(self, context: Context):
		if showOnlyThis('Model_Face'):
			return {"FINISHED"}
		else:
			self.report({'ERROR'}, "No Face Data Found")
			return {"CANCELLED"}


class updateModel(Operator):
	bl_idname = 'cgal.update_cgal_model'
	bl_label = 'Update Model'
	bl_options = {"REGISTER", "UNDO"}

	def execute(self, context: Context):
		executable = context.scene.executable
		inputFileName = context.scene.inputFile
		outputFile = path.join(path.dirname(path.abspath(executable)), "output.txt")

		if bpy.context.active_object != None and bpy.context.active_object.mode == 'EDIT':
			bpy.ops.object.mode_set(mode='OBJECT')

		process = subprocess.run([executable, inputFileName, outputFile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout = process.stdout.decode("utf-8").replace('\r', '')
		stderr = process.stderr.decode("utf-8").replace('\r', '')
		print(stdout)

		if process.returncode == 0:
			self.report({'INFO'}, "External Program Completed")
			return loadModelFromFile(outputFile, context)
		else:
			self.report({'ERROR'}, stdout + "\n" + stderr)
			return {"CANCELLED"}

class loadModel(Operator):
	bl_idname = 'cgal.load_cgal_model'
	bl_label = 'Load Model'
	bl_options = {"REGISTER", "UNDO"}

	def execute(self, context: Context):
		outputFile = context.scene.outputFile

		if bpy.context.active_object != None and bpy.context.active_object.mode == 'EDIT':
			bpy.ops.object.mode_set(mode='OBJECT')

		return loadModelFromFile(outputFile, context)

# class CGALIntegration(bpy.types.Operator):
# 	bl_idname = 'mesh.cgal_integration'
# 	bl_label = 'Add Cube'
# 	bl_options = {"REGISTER", "UNDO"}

# 	def execute(self, context):
# 		bpy.ops.mesh.primitive_cube_add()
# 		return {"FINISHED"}


def register():
	bpy.utils.register_class(vertexModel)
	bpy.utils.register_class(edgeModel)
	bpy.utils.register_class(faceModel)
	bpy.utils.register_class(updateModel)
	bpy.utils.register_class(loadModel)
	bpy.utils.register_class(vertexSelector)
	bpy.utils.register_class(edgeSelector)
	bpy.utils.register_class(faceSelector)
	bpy.utils.register_class(cgalPanel)
	bpy.types.Scene.inputFile = bpy.props.StringProperty(name="InputFile", subtype='FILE_PATH')
	bpy.types.Scene.executable = bpy.props.StringProperty(name="Executable", subtype='FILE_PATH')
	bpy.types.Scene.outputFile = bpy.props.StringProperty(name="OutputFile", subtype='FILE_PATH')
	bpy.types.Scene.vertsToSelect = bpy.props.StringProperty(name="Vertices", subtype='BYTE_STRING')
	bpy.types.Scene.edgesToSelect = bpy.props.StringProperty(name="Edges", subtype='BYTE_STRING')
	bpy.types.Scene.facesToSelect = bpy.props.StringProperty(name="Faces", subtype='BYTE_STRING')

	# handle the keymap
	wm = bpy.context.window_manager

	km = wm.keyconfigs.addon.keymaps.new(name='Object Mode', space_type='EMPTY')
	kmi = km.keymap_items.new(vertexModel.bl_idname, 'P', 'PRESS', ctrl=False, shift=False)
	addon_keymaps.append((km, kmi))
	kmi = km.keymap_items.new(edgeModel.bl_idname, 'E', 'PRESS', ctrl=False, shift=False)
	addon_keymaps.append((km, kmi))
	kmi = km.keymap_items.new(faceModel.bl_idname, 'F', 'PRESS', ctrl=False, shift=False)
	addon_keymaps.append((km, kmi))


def unregister():
	bpy.utils.unregister_class(vertexModel)
	bpy.utils.unregister_class(edgeModel)
	bpy.utils.unregister_class(faceModel)
	bpy.utils.unregister_class(updateModel)
	bpy.utils.unregister_class(loadModel)
	bpy.utils.unregister_class(vertexSelector)
	bpy.utils.unregister_class(edgeSelector)
	bpy.utils.unregister_class(faceSelector)
	bpy.utils.unregister_class(cgalPanel)
	del bpy.types.Scene.inputFile
	del bpy.types.Scene.executable
	del bpy.types.Scene.outputFile
	del bpy.types.Scene.vertsToSelect
	del bpy.types.Scene.edgesToSelect
	del bpy.types.Scene.facesToSelect

	for km, kmi in addon_keymaps:
		km.keymap_items.remove(kmi)
	addon_keymaps.clear()


class cgalPanel(Panel):
	bl_idname = "panel.cgal"
	bl_label = "CGAL"
	bl_space_type = "VIEW_3D"
	bl_region_type = "TOOLS"
	# bl_category = "Tools"
	bl_category = "CGAL"

	def draw(self, context: Context):
		layout = self.layout
		row = layout.row()
		row.operator("cgal.show_cgal_vertex_model", icon="VERTEXSEL")
		row = layout.row()
		row.operator("cgal.show_cgal_edge_model", icon="EDGESEL")
		row = layout.row()
		row.operator("cgal.show_cgal_face_model", icon="FACESEL")
		col = layout.column()
		col.prop(context.scene, 'inputFile')
		col = layout.column()
		col.prop(context.scene, 'executable')
		row = layout.row()
		row.operator("cgal.update_cgal_model", icon="FILE_REFRESH")

		col = layout.column()
		col.prop(context.scene, 'outputFile')
		row = layout.row()
		row.operator("cgal.load_cgal_model", icon="FILE_REFRESH")
		

		row = layout.row()
		row.prop(context.scene, 'vertsToSelect')
		row.operator("cgal.select_given_vertices", icon="VERTEXSEL")

		row = layout.row()
		row.prop(context.scene, 'edgesToSelect')
		row.operator("cgal.select_given_edges", icon="EDGESEL")

		row = layout.row()
		row.prop(context.scene, 'facesToSelect')
		row.operator("cgal.select_given_faces", icon="FACESEL")


if __name__ == "__main__":
	register()
