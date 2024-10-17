from ovito.io import import_file, export_file
from ovito.modifiers import *
import numpy as np
import functools

def modify_pos(frame, data, pos1):
	pos2 = data.particles_.positions_

	seg_atom_nums     = [21]*24
	seg_atom_nums[0]  = 22
	seg_atom_nums[-1] = 23
	seg_atom_nums = seg_atom_nums*36
	seg_atom_nums = seg_atom_nums + [3]*27000

	for i in range(len(pos2)):
		index0 = int(np.sum(seg_atom_nums[:i]))
		index1 = int(np.sum(seg_atom_nums[:i+1]))
		pos2[i] = np.average(pos1[index0:index1,:], axis=0)

	pos2 -= np.average(pos2, axis=0)

cnc_model_file = "cnc1+water27000.pdb"
cg_model_file  = "c864+w27000.pdb"
out_name = "cg_cnc_wat.data"

pipeline1 = import_file(cnc_model_file)
pipeline2 = import_file(cg_model_file)

data1 = pipeline1.compute()
pos1 = data1.particles['Position']

modify_temp = functools.partial(modify_pos, pos1=pos1)
pipeline2.modifiers.append(modify_temp)

export_file(pipeline2, out_name, "lammps/data", atom_style="full", ignore_identifiers=True)
