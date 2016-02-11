# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>

# TODO:

bl_info = {
    'name': 'Export BF3D Format (.bf3d)',
    'author': 'Tarcontar',
    'version': (0, 0, 1),
    "blender": (2, 6, 5),
    "api": 36079,
    'location': 'File > Export > ByteForge3D (.bf3d)',
    'description': 'Export to the ByteForge3D-Format (.bf3d)',
    'warning': 'Still in Progress',
	'tracker_url': '',
    'category': 'Import'}

# To support reload properly, try to access a package var, if it's there,
# reload everything
if "bpy" in locals():
    import imp
    if 'export_bf3d' in locals():
        imp.reload(export_bf3d)
        imp.reload(struct_bf3d)

import time
import datetime
import bpy
from bpy.props import StringProperty, BoolProperty, EnumProperty
from bpy_extras.io_utils import ImportHelper, ExportHelper
		
class ExportBF3D(bpy.types.Operator, ExportHelper):
    '''Export to bf3d file format (.bf3d)'''
    bl_idname = 'export_mesh.bf3d'
    bl_label = 'Export BF3D'
    bl_options = {'UNDO'}
	
    filename_ext = '.bf3d'
    filter_glob = StringProperty(default='*.bf3d', options={'HIDDEN'})
	
    EXPORT_MODE = EnumProperty(
            name="Export Mode",
            items=(('M', "Model", "this will export all the meshes of the scene, without skeletons or animation"), 
			('H', "Hierarchy", "this will export the hierarchy tree without any geometry or animation data"), 
			('A', "Animation", "this will export the animation without any geometry data or skeletons"), 
			),
			default='M',)
		
    def execute(self, context):
        from . import export_bf3d
        keywords = self.as_keywords(ignore=("filter_glob", "check_existing", "filepath"))		

        print('Exporting file', self.filepath)
        t = time.mktime(datetime.datetime.now().timetuple())
        export_bf3d.MainExport(self.filepath, context, self, **keywords)
        t = time.mktime(datetime.datetime.now().timetuple()) - t
        print('Finished exporting in', t, 'seconds')
        return {'FINISHED'}	

		
def menu_func_export(self, context):
    self.layout.operator(ExportBF3D.bl_idname, text='ByteForge 3D (.bf3d)')

def register():
    bpy.utils.register_module(__name__)
    bpy.types.INFO_MT_file_export.append(menu_func_export)

def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.INFO_MT_file_export.remove(menu_func_export)

if __name__ == "__main__":
    register()
