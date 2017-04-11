import bpy
from sys import argv
from os import remove

if '--' in argv:
	print(argv.index('--'))
	print(argv)

	arguments = argv[argv.index('--') + 1:]
	print(arguments)
	bpy.ops.import_mesh.off(filepath=arguments[0])
	remove(arguments[0])
	# bpy.context.scene.objects.active = bpy.data.objects[arguments[0]]
	bpy.ops.object.editmode_toggle()
	bpy.ops.mesh.select_all(action='SELECT')
	# bpy.ops.mesh.quads_convert_to_tris()
	bpy.ops.mesh.normals_make_consistent(inside=False)
	bpy.ops.mesh.select_all(action='DESELECT')
	bpy.ops.object.editmode_toggle()
	bpy.ops.export_mesh.off(filepath=arguments[0])	