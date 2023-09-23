"""
将Blender生成的.obj文件转换成Tetgen的文件格式(.node, .ele, .face, .edge, 都以0为起始索引)
"""
import numpy as np
import sys
import os
def extract_node_ele(obj_file):
    with open(obj_file) as f:
        lines = f.readlines()
        node_idx, ele_idx = 0, 0
        node_list, ele_list = [], []
        for line in lines:
            if line.startswith('v'):
                content = line.split()[1:]
                pos = [float(i) for i in content]
                node_list.append([node_idx] + pos)
                node_idx += 1
            elif line.startswith('f'):
                content = line.split()[1:]
                ele = [int(i.split('/')[0]) for i in content]
                ele_list.append([ele_idx] + ele)
                ele_idx += 1
    return node_list, ele_list

def extract_surface(tetraedrons):
    surface = set()
    for tet in tetraedrons:
        for face in (  (tet[0], tet[1], tet[2]), 
                    (tet[0], tet[2], tet[3]), 
                    (tet[0], tet[3], tet[2]),
                    (tet[1], tet[3], tet[2]) ):
            if face in surface:    
                surface.remove(face)
            else:                   
                surface.add((face[0], face[1], face[2]))
    return surface


def extrac_edge(tetraedrons):
    edges = set()
    for tet in tetraedrons:
        for edge in ( (tet[0], tet[1]), (tet[0], tet[2]), (tet[0], tet[3]),(tet[1], tet[2]), (tet[1], tet[3]), (tet[2], tet[3]) ):
            edges.add((edge[0], edge[1]))
    return edges

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python convert.py <obj_file>")
        sys.exit(1)
    input_file = sys.argv[1]
    # input_file = "./liver.tmsh"
    
    node_list, ele_list = extract_node_ele(input_file)
    tetraedrons = np.array([row[1:] for row in ele_list], dtype=np.int32)
    tetraedrons.sort(axis=1)
    surfaces = extract_surface(tetraedrons)
    edges = extrac_edge(tetraedrons)
    
    file_name = os.path.basename(input_file).split('.')[0]
    
    output_dir = os.path.dirname(input_file) + "/" 
    node_file =  output_dir + file_name + '.node'
    ele_file = output_dir + file_name + '.ele'
    surface_file = output_dir + file_name + '.face'
    edges_file = output_dir + file_name + '.edge'
    
    with open(node_file, 'w') as f:
        f.write("{} 3 0 0\n".format(len(node_list)))
        for node in node_list:
            f.write("{} {} {} {}\n".format(node[0], node[1], node[2], node[3]))
            
    with open(ele_file, 'w') as f:
        f.write("{} 4 0\n".format(len(ele_list)))
        for ele in ele_list:
            f.write("{} {} {} {} {}\n".format(ele[0], ele[1], ele[2], ele[3], ele[4]))
            
    with open(surface_file, 'w') as f:
        f.write("{} 3 0\n".format(len(surfaces)))
        for i, surface in enumerate(surfaces):
            f.write("{} {} {} {}\n".format(i, surface[0], surface[1], surface[2]))
            
    with open(edges_file, "w") as f:
        f.write("{} 2\n".format(len(edges)))
        for i, edge in enumerate(edges):
            f.write("{} {} {}\n".format(i, edge[0], edge[1]))
    