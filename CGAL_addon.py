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


class vertexModel(Operator):
	bl_idname = 'mesh.show_cgal_vertex_model'
	bl_label = 'Show Vertex'
	bl_options = {"REGISTER", "UNDO"}

	def execute(self, context: Context):
		if showOnlyThis('Model_Vertex'):
			return {"FINISHED"}
		else:
			self.report({'ERROR'}, "No Vertex Data Found")
			return {"CANCELLED"}


class edgeModel(Operator):
	bl_idname = 'mesh.show_cgal_edge_model'
	bl_label = 'Show Edge'
	bl_options = {"REGISTER", "UNDO"}

	def execute(self, context: Context):
		if showOnlyThis('Model_Edge'):
			return {"FINISHED"}
		else:
			self.report({'ERROR'}, "No Edge Data Found")
			return {"CANCELLED"}


class faceModel(Operator):
	bl_idname = 'mesh.show_cgal_face_model'
	bl_label = 'Show Face'
	bl_options = {"REGISTER", "UNDO"}

	def execute(self, context: Context):
		if showOnlyThis('Model_Face'):
			return {"FINISHED"}
		else:
			self.report({'ERROR'}, "No Face Data Found")
			return {"CANCELLED"}


class updateModel(Operator):
	bl_idname = 'mesh.update_cgal_model'
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
			scene = context.scene

			for obj in bpy.data.objects:
				obj.hide = False
				obj.select = False

			for obj in bpy.data.objects:
				obj.select = True
				bpy.ops.object.delete(use_global=False)

			edgeMesh = bmesh.new()  # type: bmesh.BMesh
			faceMesh = bmesh.new()  # type: bmesh.BMesh
			verticesInVertexMesh = []  # type: List[bmesh.BMVert]

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

				for edge in edges:
					meshEdge = faceMesh.edges.get((vertexHandle[edge[0]], vertexHandle[edge[1]]))
					if meshEdge == None:
						meshEdge = edgeMesh.edges.new((vertexHandle[edge[0]], vertexHandle[edge[1]]))
					meshEdge.seam = True

				faceMesh.to_mesh(me)
				faceMesh.free()
				obj = bpy.data.objects.new("Model_Face", me)
				scene.objects.link(obj)
				obj.data.show_edge_seams = True
				print("Face Model Generated")

			return {"FINISHED"}
		else:
			self.report({'ERROR'}, stdout + "\n" + stderr)
			return {"CANCELLED"}


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
	bpy.utils.register_class(cgalPanel)
	bpy.types.Scene.inputFile = bpy.props.StringProperty(name="InputFile", subtype='FILE_PATH')
	bpy.types.Scene.executable = bpy.props.StringProperty(name="Executable", subtype='FILE_PATH')

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
	bpy.utils.unregister_class(cgalPanel)
	del bpy.types.Scene.inputFile
	del bpy.types.Scene.executable

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
		row.operator("mesh.show_cgal_vertex_model", icon="VERTEXSEL")
		row = layout.row()
		row.operator("mesh.show_cgal_edge_model", icon="EDGESEL")
		row = layout.row()
		row.operator("mesh.show_cgal_face_model", icon="FACESEL")
		col = layout.column()
		col.prop(context.scene, 'inputFile')
		col = layout.column()
		col.prop(context.scene, 'executable')
		row = layout.row()
		row.operator("mesh.update_cgal_model", icon="FILE_REFRESH")


if __name__ == "__main__":
	register()
