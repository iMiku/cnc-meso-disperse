from ovito.io import import_file, export_file
from ovito.modifiers import *
import numpy as np
#from scipy.spatial.transform import Rotation as R
import functools

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix

def generate_random_rotation_mat():
    vec1 = np.array([0, 0, 1])
    vec2 = np.random.random(3)
    return rotation_matrix_from_vectors(vec1, vec2)

def generate_linear_cncs(volume_frac=0.01, aspect_ratio=10, cnc_grid_num=3):
    box_len = (cnc_grid_num**3*aspect_ratio/volume_frac)**(1/3)
    grid_size = box_len/cnc_grid_num
    center_points = np.zeros((cnc_grid_num**3, 3))
    center_index = 0
    for i in range(cnc_grid_num):
        for j in range(cnc_grid_num):
            for k in range(cnc_grid_num):
                center_points[center_index, 0] = grid_size*(0.5 + i)
                center_points[center_index, 1] = grid_size*(0.5 + j)
                center_points[center_index, 2] = grid_size*(0.5 + k)
                center_index += 1
    bead_points = np.zeros((aspect_ratio*cnc_grid_num**3, 3))
    for i in range(center_index):
        coords = np.zeros((int(aspect_ratio), 3))
        for j in range(int(aspect_ratio)):
            coords[j, 2] = -1.0*int(aspect_ratio)/2 + 0.5 + j
        coords = np.matmul(coords, generate_random_rotation_mat())
        coords += center_points[i]
        bead_points[i*int(aspect_ratio):(i+1)*int(aspect_ratio), :] = coords

    return bead_points

def generate_angles_info(aspect_ratio=10, cnc_grid_num=3):
    chain_num = cnc_grid_num**3
    chain_len = aspect_ratio
    lines = ""
    if(chain_len < 3):
        return lines
    lines = "\n Angles\n\n"
    for i in range(chain_num):
        for j in range(chain_len-2):
            index = j + i*(chain_len-2) + 1
            center_atom_id = j + i*chain_len + 2
            line = "%d 1 %d %d %d\n"%(index, center_atom_id-1, center_atom_id, center_atom_id+1)
            lines += line
    return lines

def read_data_add_angles(file_name, out_name, aspect_ratio=10, cnc_grid_num=3):
    if(aspect_ratio<3):
        return
    anlge_num = (aspect_ratio-2)*cnc_grid_num**3
    with open(file_name, "r") as f:
        lines = f.readlines()
    with open(out_name, "w") as f:
        for line in lines:
            f.write(line)
            if(line.rstrip().endswith("bonds")):
                f.write("%d angles\n"%(anlge_num))
            if(line.rstrip().endswith("bond types")):
                f.write("1 angle types\n")
        f.write(generate_angles_info(aspect_ratio=aspect_ratio, cnc_grid_num=cnc_grid_num))


def ovito_modify_add_cncs(frame, data, coords, aspect_ratio):
    mol_id = data.particles['Molecule Identifier'][-1] + 1
    new_bond_pairs = []
    for i in range(len(coords)):
        data.particles_.add_particle(coords[i])
        data.particles_.identifiers_[-1] = data.particles.count
        data.particles_.particle_types_[-1] = 1
        data.particles_['Molecule Identifier'][-1] = mol_id
        if(i%aspect_ratio>0):
            new_bond_pairs.append( (data.particles.count-2, data.particles.count-1) )
    bonds = data.particles_.create_bonds( count=len(new_bond_pairs) )
    bonds.create_property('Topology', data=new_bond_pairs)
    #bonds.bond_types_[:] = 1

def ovito_generate_cncs(scale_factor=4.0, out_file="out.data", volume_frac=0.01, aspect_ratio=10, cnc_grid_num=3):
    coords = generate_linear_cncs(volume_frac=volume_frac, aspect_ratio=aspect_ratio, cnc_grid_num=cnc_grid_num)
    pipeline = import_file("base.data")
    modifier1 = functools.partial(ovito_modify_add_cncs, coords=coords, aspect_ratio=aspect_ratio)
    pipeline.modifiers.append(modifier1)
    pipeline.modifiers.append(ExpressionSelectionModifier(expression = 'ParticleIdentifier==1'))
    pipeline.modifiers.append(DeleteSelectedModifier())
    scale1 = AffineTransformationModifier(
          operate_on = {'particles'}, # Transform particles but not the box.
          transformation = [[scale_factor,  0,  0, 0],
                            [ 0, scale_factor,  0, 0],
                            [ 0,  0, scale_factor, 0]],
          only_selected = False)
    pipeline.modifiers.append(scale1) 
    data = pipeline.compute()
    box_len = (cnc_grid_num**3*aspect_ratio/volume_frac)**(1/3) * scale_factor
    #print(box_len/data.cell[0,0])
    print(box_len)
    scale2 = AffineTransformationModifier(
          operate_on = {'cell'}, # Transform the box.
          transformation = [[box_len/data.cell[0,0]+1e-12,  0,  0, 0],
                            [ 0, box_len/data.cell[1,1]+1e-12,  0, 0],
                            [ 0,  0, box_len/data.cell[2,2]+1e-12, 0]],
          only_selected = False)
    pipeline.modifiers.append(scale2) 
    export_file(pipeline, out_file, "lammps/data", atom_style="full", ignore_identifiers=True)

if __name__ == '__main__':
    volume_fracs = [0.0025*x for x in range(1,7)]
    aspect_ratios= [(20*x + 10) for x in range(0,5)]
    grid_num = 6
    for vol_frac in volume_fracs:
        for aspect_ratio in aspect_ratios:
            out_name = "gn_%d_vf_%.04f_ar_%03d.data"%(grid_num, vol_frac, int(aspect_ratio) )
            ovito_generate_cncs(volume_frac=vol_frac, aspect_ratio=aspect_ratio, cnc_grid_num=grid_num, out_file=out_name)
            read_data_add_angles(file_name=out_name,out_name=out_name, aspect_ratio=aspect_ratio, cnc_grid_num=grid_num)

