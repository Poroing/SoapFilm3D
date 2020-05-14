import sys
import argparse
import pathlib

def setMat_bubble(mesh, C1, roughness, transparency=0.1, contrast=-0.5):
    mat = bpy.data.materials.new('MeshMaterial')
    mesh.data.materials.append(mat)
    mesh.active_material = mat
    mat.use_nodes = True
    tree = mat.node_tree

    # init color node
    C1 = initColorNode(tree, C1)
    C1.inputs['Contrast'].default_value = contrast 
    print(C1.inputs.keys())

    # construct car paint node
    LW = tree.nodes.new('ShaderNodeLayerWeight')
    LW.inputs['Blend'].default_value = transparency
    GLO = tree.nodes.new('ShaderNodeBsdfGlossy')
    GLO.inputs['Roughness'].default_value = roughness
    TRAN = tree.nodes.new('ShaderNodeBsdfTransparent')
    tree.links.new(C1.outputs['Color'], TRAN.inputs['Color'])
    tree.links.new(C1.outputs['Color'], GLO.inputs['Color'])

    MIX = tree.nodes.new('ShaderNodeMixShader')
    tree.links.new(LW.outputs['Facing'], MIX.inputs['Fac'])
    tree.links.new(TRAN.outputs['BSDF'], MIX.inputs[1])
    tree.links.new(GLO.outputs['BSDF'], MIX.inputs[2])

    tree.links.new(MIX.outputs[0], tree.nodes['Material Output'].inputs['Surface'])

def setEnvironmentTexture(texture_file_name):
    import bpy_extras
    world_tree = bpy.context.scene.world.node_tree
    environment_texture_node = world_tree.nodes.new('ShaderNodeTexEnvironment')
    environment_image = bpy_extras.image_utils.load_image(texture_file_name)
    environment_texture_node.image = environment_image
    world_tree.links.new(
            environment_texture_node.outputs['Color'],
            world_tree.nodes['World Output'].inputs['Surface'])

argument_parser = argparse.ArgumentParser()
argument_parser.add_argument('output_image')
argument_parser.add_argument('input_obj')
argument_parser.add_argument('--image-size-x', type=int, default=1000)
argument_parser.add_argument('--image-size-y', type=int, default=1000)
argument_parser.add_argument('--samples', type=int, default=50)
argument_parser.add_argument('--scale', type=float, default=1)
argument_parser.add_argument('--mesh-translation-x', type=float, default=0)
argument_parser.add_argument('--mesh-translation-y', type=float, default=0)
argument_parser.add_argument('--mesh-translation-z', type=float, default=0)
argument_parser.add_argument('--number-subdivisions', type=int, default=0)
argument_parser.add_argument('--environment_texture')
argument_parser.add_argument('--path-to-blender-toolbox', required=True)
args = argument_parser.parse_args(sys.argv[sys.argv.index('--') + 1:])

sys.path.append(args.path_to_blender_toolbox)
from include import *
print(sys.path)

## initialize blender
imgRes_x = args.image_size_x
imgRes_y = args.image_size_y
numSamples = args.samples
exposure = 1.5
blenderInit(imgRes_x, imgRes_y, numSamples, exposure)
bpy.ops.object.shade_smooth()

location = (args.mesh_translation_x, args.mesh_translation_y, args.mesh_translation_z)
rotation = (0, 0, 0)
scale = (args.scale,args.scale,args.scale)
mesh = readOBJ(args.input_obj, location, rotation, scale)

subdivision_surface_modifier = mesh.modifiers.new('Subdivider', 'SUBSURF')
subdivision_surface_modifier.render_levels = args.number_subdivisions
subdivision_surface_modifier.quality = 5

RGBA = (1, 1, 1, 1)
meshColor = colorObj(RGBA, 0.5, 1.0, 1.0, 0.0, 2.0)
setMat_bubble(mesh, meshColor, 0.1)

camLocation = (2,2,2)
lookAtLocation = (0,0,0.5)
focalLength = 45
cam = setCamera(camLocation, lookAtLocation, focalLength)
bpy.context.object.data.type = 'ORTHO'

if args.environment_texture is not None:
    setEnvironmentTexture(args.environment_texture)

## set light
## Option1: Three Point Light System (recommended)
#setLight_threePoints(radius=4, height=10, intensity=1700, softness=6, keyLoc='left')
## Option2: simple sun light
lightAngle = (-15,-34,-155) 
strength = 2
shadowSoftness = 0.1
sun = setLight_sun(lightAngle, strength, shadowSoftness)

## set ambient light
setLight_ambient(color=(0.1,0.1,0.1,1)) # (UI: Scene > World > Surface > Color)

## save blender file so that you can adjust parameters in the UI
#bpy.ops.wm.save_mainfile(filepath='./test.blend')

## save rendering
print(f'Rendering: {args.output_image}')
renderImage(args.output_image, cam)
print('Stoped rendering')

