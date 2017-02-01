import bpy
import bmesh

VERTS = 1
EDGES = 2
FACES = 3


def addShapeKey(obj, frame: int, mode: int):
	kn = "stage{0:0>5}".format(frame)
	sk = obj.shape_key_add(kn)  #, from_mix=False)
	bm = bmesh.new()
	bm.from_mesh(obj.data)
	bm.verts.ensure_lookup_table()
	sl = bm.verts.layers.shape.get(kn)
	for i in range(frame, frame + 1):
		if i % 100 == 0:
			print(i)
		for u in range(mode):
			bm.verts[mode * i + u][sl] *= 0
	bm.to_mesh(obj.data)


obj = bpy.context.active_object
scene = bpy.context.scene
bpy.ops.object.shape_key_remove(all=True)
sk0 = obj.shape_key_add("Basis")
print(sk0)
sk0 = obj.data.shape_keys
print(sk0)
sk0.use_relative = False

frameScale = 1

currentFrame = 1
obj.data.shape_keys.eval_time = 10
obj.data.shape_keys.keyframe_insert('eval_time', frame=currentFrame)
currentFrame = currentFrame + frameScale

for i in range(len(obj.data.vertices)):
	addShapeKey(obj, i, FACES)
	obj.data.shape_keys.eval_time = 10 * (i + 2)
	obj.data.shape_keys.keyframe_insert('eval_time', frame=currentFrame)
	currentFrame = currentFrame + frameScale

scene.frame_end = currentFrame
bpy.ops.object.shape_key_retime()
