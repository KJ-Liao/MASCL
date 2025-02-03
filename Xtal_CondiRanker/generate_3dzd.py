# Copyright (C) 2021 Tunde Aderinwale, Daisuke Kihara, and Purdue University
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# convert pdb file to triangulation data and do 3dzd calculation

import os
import sys
import glob
import numpy as np
from multiprocessing import Pool

def plytoobj(filename):
	obj_filename = filename[:-4] + '.obj'
	obj_file = open(obj_filename, 'w')

	with open(filename) as ply_file:
		ply_file_content = ply_file.read().split('\n')[:-1]

		for content in ply_file_content:
			content_info = content.split()
			if len(content_info) == 6:
				vertex_info = 'v ' + ' '.join(content_info[0:3])
				obj_file.write(vertex_info + '\n')
			elif len(content_info) == 7:
				vertex1, vertex2, vertex3 = map(int, content_info[1:4])
				vertex1, vertex2, vertex3 = vertex1 + 1, vertex2 + 1, vertex3 + 1
				face_info = 'f ' + str(vertex1) + ' ' + str(vertex2) + ' ' + str(vertex3)
				obj_file.write(face_info + '\n')

		obj_file.close()

def pdbtoinv(filename,output_dir):
	ply_dir = output_dir
	fileid = filename.split('/')[-1].split('.')[0]

	surf_command = './bin/EDTSurf -i ' + filename + ' -h 2 -o ' + ply_dir + fileid
	os.system(surf_command)

	plyfile = ply_dir + fileid + '.ply'
	if os.path.isfile(plyfile):
		plytoobj(plyfile)

		# generate 3dzd
		obj_file = ply_dir + fileid + '.obj'
		cp_command = 'cp ' + obj_file + ' ./' + fileid + '.obj'
		os.system(cp_command)

		grid_command = './bin/obj2grid -g 64  ./' + fileid + '.obj'
		os.system(grid_command)

		inv_command = './bin/map2zernike ' + fileid + '.obj.grid -c 0.5 --save-moments'
		os.system(inv_command)
		
		mv_command = 'mv ' + fileid + '.* ' + ply_dir
		os.system(mv_command)

pdb_file = sys.argv[1]
output_dir = sys.argv[2]
pdbtoinv(pdb_file,output_dir)
