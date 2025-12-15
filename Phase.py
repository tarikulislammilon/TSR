
# Program to calculate shape based structural similarities between two molecules using phase method

import itertools
import numpy as np
import math
import Bio.PDB
from Bio.PDB import PDBParser
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from sklearn.metrics.pairwise import cosine_similarity
from itertools import combinations
import pandas as pd
#from rdkit import Chem
import warnings
warnings.filterwarnings('ignore')
# Get the global periodic table object
#pt = Chem.GetPeriodicTable()

# Get the vdW radius for a specific element (e.g., Carbon)
#carbon_vdw_radius = pt.GetRvdw(8)
#print(carbon_vdw_radius)

atom_radius_dict={'C':1.7,'O':1.55,'N':1.6,'H':1.2}

df=pd.read_csv('sample_detail_tsr_11_seed10_size20 copy.csv')
PDB_list=df['protein'].to_list()
Chain=df['chain'].to_list()
drug_name=df['drug_name'].to_list()
drug_id=df['drug_id'].to_list()
Group=df['group'].to_list()

output_file=open('similarity_drugs.txt','w')

pairs=itertools.combinations(PDB_list,2)

def Distance_matrix(arr):
    diff = arr[:, np.newaxis, :] - arr[np.newaxis, :, :]
    squared_diff_sum = np.sum(diff ** 2, axis=-1)
    distance_matrix = np.sqrt(squared_diff_sum)
    return distance_matrix

def Cosine_Similarity(arr1,arr2):
    similarity_matrix = cosine_similarity(arr1, arr2)
    return similarity_matrix

def tri_alignment(Vertices_1,Vertices_2,structure1,structure2):
    # The centriod of the triangles is at origin
    # Accepts coordinates of the vertices of triangles

    Cen_triad1, Cen_triad2 = np.mean(Vertices_1, axis=0), np.mean(Vertices_2, axis=0)

    # Translates the triad at its centroid
    triad1_transformation, triad2_transformation = Vertices_1 - Cen_triad1, Vertices_2 - Cen_triad2

    structure2_transformation = structure2 - Cen_triad2
    structure1_transformation = structure1 - Cen_triad1

    # Rotate triangles to x_y plane
    u1, v1 = triad1_transformation[1] - triad1_transformation[0], triad1_transformation[2] - triad1_transformation[0]
    u2, v2 = triad2_transformation[1] - triad2_transformation[0], triad2_transformation[2] - triad2_transformation[0]

    W1 = np.cross(u1, v1)
    W2 = np.cross(u2, v2)

    magnitude_w1 = np.linalg.norm(W1)
    magnitude_w2 = np.linalg.norm(W2)

    unit_Vector_w1 = W1 / magnitude_w1
    unit_Vector_w2 = W2 / magnitude_w2

    z = np.array([0, 0, 1])

    theta1 = math.acos((np.dot(unit_Vector_w1, z)))
    theta2 = math.acos((np.dot(unit_Vector_w2, z)))

    axis1 = np.cross(unit_Vector_w1, z)
    axis2 = np.cross(unit_Vector_w2, z)

    axis1 = axis1 / np.linalg.norm(axis1)
    axis2 = axis2 / np.linalg.norm(axis2)

    I = np.identity(3)
    # Identity matrix
    K1 = np.array([
        [0, -axis1[2], axis1[1]],
        [axis1[2], 0, -axis1[0]],
        [-axis1[1], axis1[0], 0]
    ])  # Skew-symmetric matrix of the axis

    K2 = np.array([
        [0, -axis2[2], axis2[1]],
        [axis2[2], 0, -axis2[0]],
        [-axis2[1], axis2[0], 0]
    ])  # Skew-symmetric matrix of the axis

    Rotation_matrix1 = I + math.sin(theta1) * K1 + (1 - math.cos(theta1)) * np.dot(K1, K1)
    Rotation_matrix2 = I + math.sin(theta2) * K2 + (1 - math.cos(theta2)) * np.dot(K2, K2)

    inv_Rotation_matrix1=Rotation_matrix1.T

    Rotated_triad1 = np.dot(triad1_transformation, Rotation_matrix1.T)
    Rotated_triad2 = np.dot(triad2_transformation, Rotation_matrix2.T)

    Rotated_structure1 = np.dot(structure1_transformation, Rotation_matrix1.T)
    Rotated_structure2 = np.dot(structure2_transformation, Rotation_matrix2.T)

    #print(Rotated_structure2,Rotated_structure2)

    # Rotate a vertex of triad b along positive x axis
    angle = -(math.atan2(Rotated_triad2[0][1], Rotated_triad2[0][0]))
    sin_angle = math.sin(angle)
    cos_angle = math.cos(angle)

    rotation_matrix_vertices = np.array([

        [cos_angle, -sin_angle, 0],
        [sin_angle, cos_angle, 0],
        [0, 0, 1]
    ])

    #Rotated_triad2_Vertics_xaxis = np.dot(Rotated_triad2, rotation_matrix_vertices.T)
    Rotated_triad2_Vertics_xaxis = Rotated_triad2
    # Rotated_triad2_Vertics_xaxis=np.round(Rotated_triad2_Vertics_xaxis,3)

    Rotated_structure2_Vertics_xaxis = np.dot(Rotated_structure2, rotation_matrix_vertices.T)
    #print(f"Rotated_triad2_Vertics_xaxis:{Rotated_triad2_Vertics_xaxis}")

    r_b1 = np.linalg.norm(Rotated_triad2_Vertics_xaxis[0])
    r_b2 = np.linalg.norm(Rotated_triad2_Vertics_xaxis[1])
    r_b3 = np.linalg.norm(Rotated_triad2_Vertics_xaxis[2])


    beta_12 = np.arctan2(Rotated_triad2_Vertics_xaxis[0][1], Rotated_triad2_Vertics_xaxis[0][0])
    beta_23 = np.arctan2(Rotated_triad2_Vertics_xaxis[1][1], Rotated_triad2_Vertics_xaxis[1][0])
    beta_13 = np.arctan2(Rotated_triad2_Vertics_xaxis[2][1], Rotated_triad2_Vertics_xaxis[2][0])

    rB=np.array([r_b1,r_b2,r_b3])
    beta=np.array([beta_12,beta_23,beta_13])


    #xa3,xa2,xa1,ya3,ya2,ya1 = Rotated_triad1[0][0], Rotated_triad1[1][0],Rotated_triad1[2][0],Rotated_triad1[0][1],Rotated_triad1[1][1],Rotated_triad1[2][1]
    X = np.sum(rB * (Rotated_triad1[:, 0] * np.cos(beta) + Rotated_triad1[:, 1] * np.sin(beta)))
    Y = np.sum(rB * (Rotated_triad1[:, 1] * np.cos(beta) - Rotated_triad1[:, 0] * np.sin(beta)))


    theta = np.arctan2(Y, X)

    # Rotate A by theta
    R = np.array([[np.cos(theta), -np.sin(theta), 0],
                  [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])

    Rotated_triad2_triad1 = Rotated_triad1 @ R.T
    Rotated_structure1    = Rotated_structure1 @ R.T

    return Rotated_structure1,Rotated_structure2

    #print(f"Rotated_triad2_triad1:{Rotated_triad2_triad1}")


def atom_neigh_enviroment(atom_dist_matrix):

    freq_arr=[]
    for atom_dist in atom_dist_matrix:
        freq = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
        #print(f'atom_dist{atom_dist}')
        for dist in atom_dist:
            if dist != 0 and dist <= 7:
                point1 = int(dist)
                point2 = int(dist) + 1
                if dist % 1 != 0:
                    side_tri1 = point1 - (dist - 1)
                    side_tri2 = (dist + 1) - point2
                    area_tri1 = 0.5 * side_tri1 * side_tri1
                    area_tri2 = 0.5 * side_tri2 * side_tri2
                    if dist >= 1:
                        freq[point1 - 1] += area_tri1
                    if dist <= 6:
                        freq[point2] += area_tri2
                    freq[point1] += 1 - (area_tri1 + area_tri2)
                else:
                    freq[point1 - 1] += 0.5
                    freq[point1] += 0.5
        #print(freq)
        freq_arr.append(list(freq.values()))

    return freq_arr

def Self_volume_overlap(PBD1_Coord,PBD2_Coord,PDB1_atom_types,PDB2_atom_types):
    #print(f"PBD1_Coord:{PBD1_Coord}\n",f"PBD2_Coord:{PBD2_Coord}\n")
    #print(f"PDB1_atom_types:{PDB1_atom_types}\n",f"PDB2_atom_types:{PDB2_atom_types}")
    Sum_o_ab=0
    for i in range(len(PBD1_Coord)):
        for j in range(len(PBD2_Coord)):
            r_a = atom_radius_dict[PDB1_atom_types[i]]
            r_b = atom_radius_dict[PDB2_atom_types[j]]
            d_ab = np.linalg.norm(PBD1_Coord[i] - PBD2_Coord[j])
            if i!=j:
                if d_ab < r_a+r_b:
                    if d_ab <= abs(r_a - r_b):
                        o_ab = (4 / 3) * math.pi * (min(r_a, r_b)) ** 3
                        Sum_o_ab += o_ab
                    else:
                        o_ab = (math.pi / 12) * ((r_a + r_b - d_ab) ** 2) * (
                                d_ab + (2 * (r_a + r_b)) - ((3 / d_ab) * ((r_a - r_b) ** 2)))
                        Sum_o_ab += o_ab

            else:
                o_ab=(4/3)*math.pi*(r_a**3)
                Sum_o_ab+= o_ab
            #print(f"r_a:{r_a}_{PDB1_atom_types[i]}\n",f"r_b:{r_b}_{PDB2_atom_types[j]}\n",f"PBD1_Coord:{PBD1_Coord[i]}\n",f"PBD2_Coord:{PBD2_Coord[j]}\n",f"d_ab:{d_ab}\n",f"o_ab:{o_ab}\n",f"sum_o_ab:{Sum_o_ab}\n")
    return Sum_o_ab

def Volume_overlap(PBD1_Coord,PBD2_Coord,PDB1_atom_types,PDB2_atom_types):
    #print(f"PBD1_Coord:{PBD1_Coord}\n",f"PBD2_Coord:{PBD2_Coord}\n")
    #print(f"PDB1_atom_types:{PDB1_atom_types}\n",f"PDB2_atom_types:{PDB2_atom_types}")

    Sum_o_ab=0
    for i in range(len(PBD1_Coord)):
        for j in range(len(PBD2_Coord)):
            r_a = atom_radius_dict[PDB1_atom_types[i]]
            r_b = atom_radius_dict[PDB2_atom_types[j]]
            d_ab = np.linalg.norm(PBD1_Coord[i] - PBD2_Coord[j])
            if d_ab < r_a+r_b:
                if d_ab<=abs(r_a-r_b):
                    o_ab=(4/3)*math.pi*(min(r_a,r_b))**3
                else:
                    o_ab = (math.pi / 12) * ((r_a + r_b - d_ab) ** 2) * (
                            d_ab + (2 * (r_a + r_b)) - ((3 / d_ab) * ((r_a - r_b) ** 2)))

                Sum_o_ab += o_ab
        #print(f"r_a:{r_a}_{PDB1_atom_types[i]}\n",f"r_b:{r_b}_{PDB2_atom_types[j]}\n",f"PBD1_Coord:{PBD1_Coord[i]}\n",f"PBD2_Coord:{PBD2_Coord[j]}\n",f"d_ab:{d_ab}\n",f"o_ab:{o_ab}\n",f"sum_o_ab:{Sum_o_ab}\n")
    return Sum_o_ab

def PDB_File_Parser(Structure,Chain_,Residue_,Resisdue_ID):
    model = Structure[0]
    Coord_PBD = []
    atom_types_PBD = []
    for chain in model:
        if chain.id == Chain_:
            for residue in chain:
                if str(residue)[9:12] == Residue_:
                    numeric_filter = filter(str.isdigit, str(residue.id))
                    Res_Id = "".join(numeric_filter)  # Residue ID
                    if int(Res_Id) ==Resisdue_ID:
                        for atom in residue:
                            atom_types_PBD.append(atom.get_id()[0])
                            Coord_PBD.append(atom.get_coord())

    return Coord_PBD,atom_types_PBD

def RMSD(A,B):
    #diff = arr[:, np.newaxis, :] - arr[np.newaxis, :, :]
    diff=A[:, np.newaxis, :] - B[np.newaxis, :, :]
    squared_diff_sum = np.sum(diff ** 2, axis=-1)
    distance_matrix = np.sqrt(squared_diff_sum)
    min_dist=np.min(distance_matrix,axis=1)
    rmsd=np.sqrt(np.mean(min_dist**2))
    return rmsd


Similarity_matrix=np.zeros((len(PDB_list),len(PDB_list)))
Index_list=[]
for i in range(len(PDB_list)):
    Index_list.append(f"{PDB_list[i]}_{Chain[i]}_{drug_name[i]}_{drug_id[i]}_{Group[i]}")
    for j in range(i+1,len(PDB_list)):

        PDB1, PDB2 = PDB_list[i], PDB_list[j]
        #print(PDB1,PDB2)

        PDB_File_path1 = f"/Users/c00479477/Desktop/PycharmProjects/pythonProject/PDB_data/{PDB1}.pdb"
        PDB_File_path2 = f"/Users/c00479477/Desktop/PycharmProjects/pythonProject/PDB_data/{PDB2}.pdb"

        p = Bio.PDB.PDBParser()

        Structure1 = p.get_structure('PrimaryStructureChain', PDB_File_path1)
        Structure2 = p.get_structure('PrimaryStructureChain', PDB_File_path2)

        Coord_PBD1, atom_types_PBD1 = PDB_File_Parser(Structure1,Chain[i],drug_name[i],drug_id[i])
        Coord_PBD2, atom_types_PBD2 = PDB_File_Parser(Structure2,Chain[j],drug_name[j],drug_id[j])

        #print(type(Chain[i]),type(drug_name[i]),type(drug_id[i]))

        Coord_PBD1 = np.array(Coord_PBD1)
        Coord_PBD2 = np.array(Coord_PBD2)

        Distance_matrix_PBD1 = Distance_matrix(Coord_PBD1)
        Distance_matrix_PBD2 = Distance_matrix(Coord_PBD2)

        atom_neigh_PDB1=atom_neigh_enviroment(Distance_matrix_PBD1)
        atom_neigh_PDB2=atom_neigh_enviroment(Distance_matrix_PBD2)

        Cosine_similarity_matrix=Cosine_Similarity(np.array(atom_neigh_PDB1),np.array(atom_neigh_PDB2))
        #print(Cosine_similarity_matrix)


        max_indices_per_row = np.argmax(Cosine_similarity_matrix, axis=0)
        #print(max_indices_per_row)
        max_values_per_row = np.max(Cosine_similarity_matrix, axis=0)
        #print(max_values_per_row)
        max_values_indexed_elements = tuple(enumerate(max_indices_per_row))
        #print(max_values_indexed_elements)

        max_values_dict={}
        count=0
        for max_values in max_values_indexed_elements:
            max_values_dict[max_values]=max_values_per_row[count]
            count+=1
        #print(max_values_dict)
        sorted_max_values= dict(sorted(max_values_dict.items(), key=lambda item: item[1], reverse=True))
        sorted_max_values_list=list(sorted_max_values.keys())
        #print(len(sorted_max_values_list))

        All_overlaps_AB=[]
        for l in range(len(sorted_max_values_list)):
            for m in range(l+1,len(sorted_max_values_list)):
                for n in range(m+1,len(sorted_max_values_list)):

                    b1,b2,b3,a1,a2,a3=sorted_max_values_list[l][0],sorted_max_values_list[m][0],sorted_max_values_list[n][0],sorted_max_values_list[l][1],sorted_max_values_list[m][1],sorted_max_values_list[n][1]
                    #print(b1,b2,b3,a1,a2,a3)
                    b1_coord,b2_coord,b3_coord,a1_coord,a2_coord,a3_coord=Coord_PBD2[b1],Coord_PBD2[b2],Coord_PBD2[b3],Coord_PBD1[a1],Coord_PBD1[a2],Coord_PBD1[a3]
                    #print(b1_coord,b2_coord,b3_coord,a1_coord,a2_coord,a3_coord)
                    if b1 != b2 != b3 and a1 != a2 != a3:

                        b_12, b_23, b_13 = np.linalg.norm(b1_coord-b2_coord),   np.linalg.norm(b2_coord-b3_coord),   np.linalg.norm(b1_coord-b3_coord)
                        a_12, a_23, a_13 = np.linalg.norm(a1_coord - a2_coord), np.linalg.norm(a2_coord - a3_coord), np.linalg.norm(a1_coord - a3_coord)
                        #print(b_12, b_23, b_13,a_12, a_23, a_13)

                        if 0 <= abs(b_12 - a_12) <= 2 and 0 <= abs(b_23 - a_23) <=2 and 0 <= b_13 - a_13 <= 2 :

                            Triad_1_Vertices = np.array([a1_coord, a2_coord, a3_coord], dtype=float)
                            Triad_2_Vertices = np.array([b1_coord, b2_coord, b3_coord],dtype=float)

                            #print(Triad_1_Vertices,Triad_2_Vertices)
                            aligned_Coord_PDB1,aligned_Coord_PDB2=tri_alignment(Triad_1_Vertices,Triad_2_Vertices,Coord_PBD1,Coord_PBD2)
                            #Overlap_A_B=Volume_overlap(aligned_Coord_PDB1,aligned_Coord_PDB2,atom_types_PBD1,atom_types_PBD2)
                            #For RMSD values
                            Overlap_A_B=RMSD(aligned_Coord_PDB1,aligned_Coord_PDB2)
                            All_overlaps_AB.append(Overlap_A_B)
                        else:
                            All_overlaps_AB.append(0)
                    else:
                        All_overlaps_AB.append(0)
                            #print(np.mean(Coord_PBD1, axis=0), np.mean(aligned_structure_2, axis=0))
        if len(All_overlaps_AB)!=0:
            O_AA = Self_volume_overlap(Coord_PBD1, Coord_PBD1, atom_types_PBD1, atom_types_PBD1)
            O_BB = Self_volume_overlap(Coord_PBD2, Coord_PBD2, atom_types_PBD2, atom_types_PBD2)
            # print(O_AA,O_BB)
            # print(max(All_overlaps_AB))
            Similarity = max(All_overlaps_AB) / max([O_AA, O_BB])
            Similarity_matrix[i, j] = Similarity
            Similarity_matrix[j, i] = Similarity
            output_file.write(f"{PDB1}_{Chain[i]}_{drug_name[i]}_{drug_id[i]}_{Group[i]}\t{PDB2}_{Chain[j]}_{drug_name[j]}_{drug_id[j]}_{Group[j]}\t{Similarity}\n")
        else:
            print(PDB_list[i],Chain[i],drug_name[i],drug_id[i], PDB_list[j],Chain[j],drug_name[j],drug_id[j])
        #print(PDB1,PDB2,Similarity)

similarity_df = pd.DataFrame(Similarity_matrix, index=Index_list)

# 5. Save to CSV
similarity_df.to_csv('similarity_matrix.csv', index=True,header=False)
output_file.close()
print('completed_succesfully')

