#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import math

pdb_file_path = '/Users/ymlee/Desktop/research/02_tutorial_protein_polygon/RK718_full.pdb'

# 읽기 모드로 file open, lines 리스트에 각 줄 저장
# with 사용시 블록이 끝나면 자동으로 file close
with open(pdb_file_path,'r') as pdb_file:
    lines = pdb_file.readlines()
    if lines:
        last_line = lines[-1] # -1: 마지막 index
    else:
        last_line = None

if last_line != None:
    lines_num = len(lines)
    #print (lines_num)
    #print("check end: ",last_line)
else:
    print("empty file")

# splitted_lines의 각 원소는 각 줄을 split하여 저장한 문자열로 구성됨 (~2차원 array)
splitted_lines = []
for i in range(0,lines_num):
    splitted_lines.append(lines[i].split())
# print(splitted_lines)


# In[2]:


# Atom class
class Atom:
    # 생성자 정의
    def __init__(self, line): 
        self.num = line[1] # 인스턴스/클래스 변수는 기본적으로 public. 앞에 __ 붙이면 private
        self.name = line[2]
        self.aa = line[3]
        self.chain = line[4]
        self.residue = line[5]
        self.x = line[6]
        self.y = line[7]
        self.z = line[8]
        self.coord = [float(line[6]),float(line[7]),float(line[8])]
        self.symbol_charge = line[11]

    # printing elements
    def print_num(self):
        print(self.num)
    def print_name(self):
        print(self.name)
    def print_aa(self):
        print(self.aa)
    def print_chain(self):
        print(self.chain)
    def print_residue(self):
        print(self.residue)
    def print_coord(self):
        print(self.coord)
    def print_symbol_charge(self):
        print(self.symbol_charge)
    
    '''
    def print_num_nolinebreak(self):
        print(self.num, end = "\t")
    def print_chain_nolinebreak(self):
        print(self.chain, end = "\t")
    def print_residue_nolinebreak(self):
        print(self.residue, end = "\t")
    '''
    def print_num_chain_residue_nolinebreak(self):
        print(self.num,"\t",self.chain,"\t",self.residue,end = "\t")


# In[3]:


# Global variables
atom_list = [] # List of Atom classes
selected_coords = [] # 원하는(출력한) coordinates를 하나로 합쳐 저장해 둘 공간(list)
selected_range = [[0,0]] # 번호를 저장해서 coordinates 이외의 정보는 atom_list로부터 가져옴


# In[4]:


# Initializing(appending) atom_list 
count_TER = -1 # TER 개수를 셈, 마지막에서 세 번째 line의 chain 확인으로 세는 것도 가능

for j in range(0,lines_num):
    if splitted_lines[j][0] == 'ATOM':
        atom_list.append(Atom(splitted_lines[j]))
    else:
        count_TER += 1
#print(atom_list)


# In[5]:


# Functions
# Initializing(clearing) selected_coords
def init_selected_coords(selected_coords): # (!!!)함수 내에서 호출되는 경우 global variable이더라도 매개변수로 들어가야 함
    if len(selected_coords) != 0:
        selected_coords = [] # https://velog.io/@jyunxx/problem-solvingPython-리스트의-clear-함수와-초기
        selected_range = [[0,0]]

# 전체 coord 출력
def all_coords(lines_num,count_TER):
    init_selected_coords(selected_coords) # 여러 range 추가/삭제가 필요한 경우 따로 호출
    for k in range(0,lines_num-count_TER-1):
        #atom_list[k].print_coord()
        selected_coords.append(atom_list[k].coord)
    selected_range[0] = [1,lines_num-count_TER-1] # selected_range = [[]]로 사용하면 함수 밖에서 값이 사라짐 (local var로 초기화되는 듯함)
    return selected_coords

#all_coords(lines_num,count_TER)
#print(selected_range[0])


# In[6]:


# 원자번호 범위 [start,end] 지정하면 coord 출력
def range_coords(start,end,lines_num,count_TER):
    if (start<1) or (start>lines_num-count_TER-1) or (end<1) or (end>lines_num-count_TER-1):
        print("out of range")
    elif start>end:
        print("wrong range")
    else:
        init_selected_coords(selected_coords)
        print("atom number | coordinate")
        for k in range(start-1,end):
            atom_list[k].print_num_chain_residue_nolinebreak()
            atom_list[k].print_coord()
            selected_coords.append(atom_list[k].coord)
        selected_range = [[start,end]]
    return 1
    
#range_coords(1000,2000,lines_num,count_TER)
#print(temp_coords)


# In[7]:


# chain 지정하면 coord 출력 (select X)
def chain_coords(chain_name,lines_num,count_TER):
    _coords = []
    for k in range(0,lines_num-count_TER-1):
        if atom_list[k].chain == chain_name:
            _coords.append(atom_list[k].coord)
    return _coords

#chain_coords("A",lines_num,count_TER)


# In[8]:


# chain, residue 범위 지정하면 coord 출력
def residue_coords(chain_name,lines_num,count_TER):
    last_atom_num_in_chain = 0
    for k in range(0,lines_num-count_TER-1): # python에서는 for/while 중 어떤 게 효율적일까
        if atom_list[k].chain > chain_name:
            break
        last_atom_num_in_chain += 1
    #print(last_atom_num_in_chain)
    residues_num = int(atom_list[last_atom_num_in_chain-1].residue)

    # 해당 chain 내 residue의 개수를 출력, 시작/끝 residue 번호를 입력받음
    print("The number of residues in chain ",chain_name,": ",residues_num)
    start = int(input("Enter the number of the first residue: "))
    end = int(input("Enter the number of the last residue: "))
    if (start<1) or (start>residues_num) or (end<1) or (end>residues_num):
        print("out of range")
    elif start>end:
        print("wrong range")
    else:
        init_selected_coords(selected_coords)
        temp_range = []
        print("atom number | chain name | residue number | coordinate")
        for k in range(0,lines_num-count_TER-1):
            if (atom_list[k].chain == chain_name) and (int(atom_list[k].residue) >= start) and (int(atom_list[k].residue) <= end):
                atom_list[k].print_num_chain_residue_nolinebreak()
                atom_list[k].print_coord()
                selected_coords.append(atom_list[k].coord)
                temp_range.append(int(atom_list[k].num))
        selected_range = [[temp_range[0],temp_range[len(temp_range)-1]]]
        print(selected_range)
    return 1
    
#residue_coords('A',lines_num,count_TER)


# In[9]:


# temp_coords를 np.array로 변경
# https://codingdog.tistory.com/entry/list를-numpy-array로-바꾸고-numpy-array를-list로-바꾸는-방법을-알아봅시다
def gen_nparray(temp_coords):
    nparray = np.array(temp_coords)
    #print(nparray)
    return nparray

# numpy 이용하여 coord 이동
def translation(selected_coords,x,y,z):
    coords_nparray = np.array(selected_coords)
    shift_amount = np.array(([float(x),float(y),float(z)]))
    shifted_coords_nparray = coords_nparray + shift_amount
    shifted_coords = shifted_coords_nparray.tolist() # from numpy array to list
    return shifted_coords

#print(translation(selected_coords,1000,1000,1000))
# 한 자리 수 더할 때 부동소수점 오류 발생함


# In[10]:


# 추가) selected_coords가 수정된 경우 전체 coords에 덮어쓰기 (warning 포함)
#def overwrite_selected_coords():
# 추가) 특정 residue/coord 삭제 기능 (lines_num, count_TER 수정해야 함)
# 추가) 여러 개의 chain 각각 범위 지정 가능하게


# In[11]:


# pdb 형식대로 string 이어붙여 파일 저장
# https://datascienceschool.net/01%20python/02.04%20파이썬의%20문자열%20형식화.html
def save_new_coords(coords,filename): # coords: float list type의 매개변수 (selected_coords, shifted_coords (from parallel_trans()), ...)
    alter_location_indicator = " "
    code_for_insertion_of_residues = " "
    occupancy = 1.00
    temp_factor = 0.00
    #print(selected_range[0][0])
    current_chain = "A"
    with open(filename,"w") as new_pdb_file:
        for n in range(selected_range[0][-1]):
            #print(n)
            num = atom_list[n].num
            name = atom_list[n].name
            if(name[0]<"A"):
                name1 = name[0:1]
                name2 = name[1:]
            else:
                name1 = " "
                name2 = name
            aa = atom_list[n].aa
            chain = atom_list[n].chain
            residue = atom_list[n].residue
            x = coords[n-selected_range[0][0]+1][0]
            y = coords[n-selected_range[0][0]+1][1]
            z = coords[n-selected_range[0][0]+1][2]
            symbol_charge = atom_list[n].symbol_charge
            if(len(symbol_charge)>2): # split symbol_charge into symbol and charge
                charge = symbol_charge[-2:]
                symbol = symbol_charge.replace(charge,"") # 기존 문자열은 수정 X
            else:
                charge = "  "
                symbol = symbol_charge
            
            #print("ATOM      1  N   ASP A   1      38.927  13.356  30.662  1.00  0.00           N  ") #형식 확인용
            #print("%-6s%5s %1s%-3s%1s%-3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"
            #    % ("ATOM",num,name1,name2,alter_location_indicator,aa,chain,residue,code_for_insertion_of_residues
            #,x,y,z,occupancy,temp_factor,symbol,charge))
            
            if current_chain < chain: # chain이 끝나면 TER 출력
                new_pdb_file.write("TER  \n")
                current_chain = chain
            new_pdb_file.write("%-6s%5s %1s%-3s%1s%-3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n"
                    % ("ATOM",num,name1,name2,alter_location_indicator,aa,chain,residue,code_for_insertion_of_residues
            ,x,y,z,occupancy,temp_factor,symbol,charge))
            
        new_pdb_file.write("TER  \nEND\n") 

#save_new_coords(selected_coords)


# In[12]:


### Functions from quatMath.py, quatMath_Michael.py (can be slightly modified, marked)
## rotate around axis (x, y, z) by angle theta (degree)
def getQuaternion(x, y, z, theta):
    ## Step 1. Normalization
    norm = math.sqrt(x*x + y*y + z*z)
    x /= norm
    y /= norm
    z /= norm

    ## Step 2. Generate quaternion for rotation
    q0 = math.cos(math.radians(theta*0.5))
    q1 = x * math.sin(math.radians(theta*0.5))
    q2 = y * math.sin(math.radians(theta*0.5))
    q3 = z * math.sin(math.radians(theta*0.5))

    return np.array([q0, q1, q2, q3])

## Returns the rotation of a vector by a given quaternion.
# input quat and vec are arrays of four numbers: [q0, q1, q2, q3]
def quatVec(quat, vec):
    v0 = (pow(quat[0],2) + pow(quat[1],2) - pow(quat[2],2) - pow(quat[3],2))*vec[0] \
            + (2.*quat[1]*quat[2]-2.*quat[0]*quat[3])*vec[1] + (2.*quat[3]*quat[1]+2.*quat[0]*quat[2])*vec[2]

    v1 = (2.*quat[1]*quat[2]+2.*quat[0]*quat[3])*vec[0] + (pow(quat[0],2)-pow(quat[1],2)+pow(quat[2],2)-pow(quat[3],2))*vec[1] \
            + (2.*quat[2]*quat[3]-2.*quat[0]*quat[1])*vec[2]

    v2 = (2.*quat[1]*quat[3]-2.*quat[0]*quat[2])*vec[0] + (2.*quat[0]*quat[1]+2.*quat[2]*quat[3])*vec[1] \
            + (pow(quat[0],2) - pow(quat[1],2) - pow(quat[2],2) + pow(quat[3],2))*vec[2]

    return np.array([v0, v1, v2])

## Rotates a list of vectors according to the quaternion quat
# input quat and vec are arrays of four numbers: [q0, q1, q2, q3]
def rotateFrame(quat, veclist):
    if (quat == [1,0,0,0]).all(): # modified (added if)
        return veclist
    else: 
        vecarr = np.asarray(veclist) # why asarray??
        new_veclist = []

        for vec in vecarr:
            newvec = quatVec(quat, vec)
            new_veclist.append(newvec)

        #new_veclist = np.asarray(new_veclist) # modified (no need)

        return new_veclist

## Returns the product of two quaternions.
# input quats are arrays of four numbers: [q0, q1, q2, q3]
def quatMultiply(quatA, quatB):
    # quaternion scalar
    qAs = quatA[0]
    # quaternion vector
    qAv = np.array([quatA[1], quatA[2], quatA[3]])
    # quaternion scalar
    qBs = quatB[0]
    # quaternion vector
    qBv = np.array([quatB[1], quatB[2], quatB[3]])

    # product scalar and vector
    qABs = qAs*qBs - np.dot(qAv, qBv)
    qABv = qAs*qBv + qBs*qAv + np.cross(qAv, qBv)

    return np.array([qABs, qABv[0], qABv[1], qABv[2]])

## Returns the angle between two vectors
def calAngle(v1,v2):
    # Angle between v1 and v2
    cosine_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.arccos(cosine_angle) # Radian
    angle = np.degrees(angle) # Degree
    return angle

## Returns the distance between two vectors
def calDist(v1,v2):
    return math.sqrt((v2[0]-v1[0])**2 + (v2[1]-v1[1])**2 + (v2[2]-v1[2])**2)

## gets a quaternion that rotates vec1 to vec2
def VecsToQuat(vec1, vec2):
    vec1=np.array(vec1)
    vec2=np.array(vec2)
    # get a unit vector that is perpendicular to vec1 and vec2
    # this will be the vector for which clockwise rotation about it, from vec1 to vec2, is between 0 and pi.
    n = np.cross(vec1, vec2)
    # if the angle between vec1 and vec2 is zero, don't do anything crazy
    if (np.sqrt(np.dot(n,n)) < 1e-12):
        return np.array([1,0,0,0])
    n = n/np.sqrt(np.dot(n,n))
    # get the angle between vec1 and vec2. the range of acos is [0,pi], so we are fine.
    theta = math.acos(np.dot(vec1, vec2)/np.sqrt(np.dot(vec1,vec1))/np.sqrt(np.dot(vec2,vec2)))
    # convert to a quaternion scalar and vector part
    qs = math.cos(theta/2.)
    qv = math.sin(theta/2.)*n
    # return as a 4D numpy array
    return np.array([qs, qv[0], qv[1], qv[2]])


# In[13]:


# Center of Mass
def CM(_coords):
    CM = [0.0,0.0,0.0]
    range_length = len(_coords)
    npcoords = gen_nparray(_coords)
    npCM = gen_nparray(CM)
    for i in range(range_length):
        npCM += npcoords[i]        
    npCM /= range_length
    CM = npCM.tolist()
    return CM
    
#CM(selected_coords)


# In[14]:


# returns arm length and arm vector(from origin to arm end 아멘)
def arm_length(_coords,arm_chain,arm_end_residue,origin=np.array([0,0,0])):
    arm_end_coords = []
    arm_end_range = []
    for atom in atom_list:
        if atom.chain == arm_chain and atom.residue == arm_end_residue:
            arm_end_range.append(int(atom.num))
    for i in arm_end_range:
        arm_end_coords.append(_coords[i-1])
    arm_end = gen_nparray(CM(arm_end_coords))
    arm_v = arm_end - origin
    arm_v_length = calDist(origin,arm_v)
    return arm_v_length, arm_v


# In[15]:


### Rotation and Translation (2D)
## 3개의 동일한 monomer로 구성된 trimer로 정육각형 만들기: 2가지 회전
def hexagon(arm_chain,arm_end_residue):
    # 임의의 팔(arm)을 지정, 벡터(arm_v) 계산
    # both origin and arm_v are numpy arrays
    # * origin을 (0,0,0)으로 옮기지 않고 기존 trimer의 CM으로 하여 배열
    origin = gen_nparray(CM(all_coords(lines_num,count_TER))) # trimer 중심 (전체 CM)
    arm_v_length, arm_v = arm_length(selected_coords,arm_chain,arm_end_residue,origin)

    # Quaternions(q1, q2) 계산
    # q1
    v1 = origin + np.array([0,-arm_v_length,0]) #arm_v_length
    q1 = VecsToQuat(arm_v,v1)
    #print(q1)
    # q2
    v2 = origin + np.array([0,arm_v_length,0])
    q2 = VecsToQuat(arm_v,v2)
    #print(q2)
    
    # 2가지 rotation 진행
    rot_coords1 = rotateFrame(q1,selected_coords)
    rot_coords2 = rotateFrame(q2,selected_coords)
    
    # 6가지 translation 후 파일들 저장 (제 1 -> 4 사분면 순서)
    l = 2*arm_v_length + 10 # unit: Angstrom, 육각형 한 변의 길이(edge)
    root3_l = math.sqrt(3)*l

    trimer1_coords = translation(rot_coords1,root3_l/2,l/2,0)
    trimer2_coords = translation(rot_coords2,root3_l/2,-l/2,0)
    trimer3_coords = translation(rot_coords1,0,-l,0)
    trimer4_coords = translation(rot_coords2,-root3_l/2,-l/2,0)
    trimer5_coords = translation(rot_coords1,-root3_l/2,l/2,0)
    trimer6_coords = translation(rot_coords2,0,l,0)
    save_new_coords(trimer1_coords,"RK718_full_hexagon_1.pdb")
    save_new_coords(trimer2_coords,"RK718_full_hexagon_2.pdb")
    save_new_coords(trimer3_coords,"RK718_full_hexagon_3.pdb")
    save_new_coords(trimer4_coords,"RK718_full_hexagon_4.pdb")
    save_new_coords(trimer5_coords,"RK718_full_hexagon_5.pdb")
    save_new_coords(trimer6_coords,"RK718_full_hexagon_6.pdb")

hexagon("C", "307")


# In[16]:


### Rotation and Translation (3D)
## 3개의 동일한 monomer로 구성된 trimer로 정육면체 만들기: (8가지 회전 -> 8가지 회전) =>  (2가지 회전 -> 2가지 회전) * 3번 대칭 (정다면체라 가능)
# CM of original trimer를 origin (0,0,0), trimer_ortho_v를 (0,0,1)로 이동하여 초기화 (return initialized coord.s)
def toOrigin(original_trimer_CM,A_CM,B_CM,C_CM):    
    init_coords = translation(selected_coords,-original_trimer_CM[0],-original_trimer_CM[1],-original_trimer_CM[2])
    #print(CM(init_coords)) # == (0,0,0)인지 확인
    init_A_coords = translation(chain_coords("A",lines_num,count_TER),-original_trimer_CM[0],-original_trimer_CM[1],-original_trimer_CM[2])
    init_B_coords = translation(chain_coords("B",lines_num,count_TER),-original_trimer_CM[0],-original_trimer_CM[1],-original_trimer_CM[2])
    init_C_coords = translation(chain_coords("C",lines_num,count_TER),-original_trimer_CM[0],-original_trimer_CM[1],-original_trimer_CM[2])
    A_CM = gen_nparray(CM(init_A_coords))
    B_CM = gen_nparray(CM(init_B_coords))
    C_CM = gen_nparray(CM(init_C_coords))
    #print(A_CM+B_CM+C_CM)
    AB_v = B_CM-A_CM
    AC_v = C_CM-A_CM
    trimer_ortho_v = np.cross(AC_v,AB_v) 
    q = np.array(VecsToQuat(trimer_ortho_v,[0,0,1])) # convex 부분이 +z 방향을 향하게
    #print(q)
    init_coords = rotateFrame(q,init_coords)
    return init_coords, A_CM, B_CM, C_CM

def cube_symm(arm_chain,arm_end_residue):
    original_trimer_CM = CM(all_coords(lines_num,count_TER))
    A_CM = np.array([0,0,0])
    B_CM = np.array([0,0,0])
    C_CM = np.array([0,0,0])
    init_coords,A_CM,B_CM,C_CM = toOrigin(original_trimer_CM,A_CM,B_CM,C_CM)
    CM_dist = calDist(A_CM,B_CM)
    arm_v_length, A_arm = arm_length(init_coords,"A",arm_end_residue) 
    arm_v_length, B_arm = arm_length(init_coords,"B",arm_end_residue) 
    arm_v_length, C_arm = arm_length(init_coords,"C",arm_end_residue) 
    arm_end_dist = calDist(A_arm,B_arm)
    # 정육면체 모서리 길이
    # l = 2*CM_dist + 10 # when using the arm CMs
    l = 2*arm_v_length + 10 # when using the arm ends
    
    # trimer_v = ((A_CM+B_CM)/2 - C_CM).tolist() # when using the arm CMs
    trimer_v = ((A_arm+B_arm)/2-C_arm).tolist() # when using the arm ends
    # 2 rotations
    # *순서상 1st, 2nd rotation이 바뀌었지만 두 회전 모두 (0,0,0)을 지나는 axis가 회전축이므로 가능
    root3 = math.sqrt(3)
    root6 = math.sqrt(6)    
    # 1st rotation (각 trimer의 arm이 향하는 방향이 정육면체의 모서리가 되도록 회전)
    # trimer_q1 = VecsToQuat(trimer_v,[-3*CM_dist/root6,-3*CM_dist/root6,0]) # when using the arm CMs
    trimer_q1 = VecsToQuat(trimer_v,[-3*arm_end_dist/root6,-3*arm_end_dist/root6,0]) # when using the arm ends
    trimer_rot1 = rotateFrame(trimer_q1,init_coords)
    # 2nd roation
    trimer1_q2 = VecsToQuat([0,0,1],[1/root3,1/root3,1/root3])
    trimer1_rot2 = rotateFrame(trimer1_q2,trimer_rot1)
    trimer2_q2 = VecsToQuat([0,0,1],[-1/root3,-1/root3,-1/root3]) 
    trimer2_rot2 = rotateFrame(trimer2_q2,trimer_rot1)

    # 2가지 trimer에 대한 translation 진행
    trimer1_1_coords = np.array(translation(trimer1_rot2,l/2,l/2,l/2))
    trimer2_1_coords = np.array(translation(trimer2_rot2,-l/2,-l/2,-l/2))
    
    # 3가지 대칭 진행 : 입체 이성질체가 되어 pymol에서 인식이 불가능한 것으로 추측
    '''
    # 부호 변경
    trimer1_2_coords = []
    trimer1_3_coords = []
    trimer1_4_coords = []
    trimer2_2_coords = []
    trimer2_3_coords = []
    trimer2_4_coords = [] 
    for i in range(selected_range[0][-1]):
        # x-y 면대칭
        trimer1_2_coords.append([trimer1_1_coords[i][0],trimer1_1_coords[i][1],-(trimer1_1_coords[i][2])])
        trimer2_2_coords.append([trimer2_1_coords[i][0],trimer2_1_coords[i][1],-(trimer2_1_coords[i][2])])
        # y-z 면대칭
        trimer1_3_coords.append([-(trimer1_1_coords[i][0]),trimer1_1_coords[i][1],trimer1_1_coords[i][2]])
        trimer2_3_coords.append([-(trimer2_1_coords[i][0]),trimer2_1_coords[i][1],trimer2_1_coords[i][2]])
        # z-x 면대칭
        trimer1_4_coords.append([trimer1_1_coords[i][0],-(trimer1_1_coords[i][1]),trimer1_1_coords[i][2]])
        trimer2_4_coords.append([trimer2_1_coords[i][0],-(trimer2_1_coords[i][1]),trimer2_1_coords[i][2]])
    '''
    # 행렬 변환
    symm_xy = np.array([[1,0,0],[0,1,0],[0,0,-1]])
    trimer1_2_coords = (trimer1_1_coords @ symm_xy.T).tolist()
    trimer2_2_coords = (trimer2_1_coords @ symm_xy.T).tolist()
    
    symm_yz = np.array([[-1,0,0],[0,1,0],[0,0,1]])
    trimer1_3_coords = (trimer1_1_coords @ symm_yz.T).tolist()
    trimer2_3_coords = (trimer2_1_coords @ symm_yz.T).tolist()
    
    symm_zx = np.array([[1,0,0],[0,-1,0],[0,0,1]])
    trimer1_4_coords = (trimer1_1_coords @ symm_zx.T).tolist()
    trimer2_4_coords = (trimer2_1_coords @ symm_zx.T).tolist()

    save_new_coords(init_coords,"RK718_full_cube_origin.pdb")
    save_new_coords(trimer1_1_coords,"RK718_full_cube_1_1.pdb")
    save_new_coords(trimer1_2_coords,"RK718_full_cube_1_2.pdb")
    save_new_coords(trimer1_3_coords,"RK718_full_cube_1_3.pdb")
    save_new_coords(trimer1_4_coords,"RK718_full_cube_1_4.pdb")
    save_new_coords(trimer2_1_coords,"RK718_full_cube_2_1.pdb")
    save_new_coords(trimer2_2_coords,"RK718_full_cube_2_2.pdb")
    save_new_coords(trimer2_3_coords,"RK718_full_cube_2_3.pdb")
    save_new_coords(trimer2_4_coords,"RK718_full_cube_2_4.pdb")

#cube_symm("C","307")


# In[19]:


def cube_symm(arm_chain,arm_end_residue):
    original_trimer_CM = CM(all_coords(lines_num,count_TER))
    A_CM = np.array([0,0,0])
    B_CM = np.array([0,0,0])
    C_CM = np.array([0,0,0])
    init_coords,A_CM,B_CM,C_CM = toOrigin(original_trimer_CM,A_CM,B_CM,C_CM)
    CM_dist = calDist(A_CM,B_CM)
    arm_v_length, A_arm = arm_length(init_coords,"A",arm_end_residue) 
    arm_v_length, B_arm = arm_length(init_coords,"B",arm_end_residue) 
    arm_v_length, C_arm = arm_length(init_coords,"C",arm_end_residue) 
    arm_end_dist = calDist(A_arm,B_arm)
    # 정육면체 모서리 길이
    # l = 2*CM_dist + 10 # when using the arm CMs
    l = 2*arm_v_length + 10 # when using the arm ends
    
    # trimer_v = ((A_CM+B_CM)/2 - C_CM).tolist() # when using the arm CMs
    trimer_v = ((A_arm+B_arm)/2-C_arm).tolist() # when using the arm ends
    # 2 rotations
    # *순서상 1st, 2nd rotation이 바뀌었지만 두 회전 모두 (0,0,0)을 지나는 axis가 회전축이므로 가능
    root3 = math.sqrt(3)
    root6 = math.sqrt(6)    
    # 1st rotation (각 trimer의 arm이 향하는 방향이 정육면체의 모서리가 되도록 회전)
    # trimer_q1 = VecsToQuat(trimer_v,[-3*CM_dist/root6,-3*CM_dist/root6,0]) # when using the arm CMs
    trimer1_13_q1 = VecsToQuat(trimer_v,[-3*arm_end_dist/root6,-3*arm_end_dist/root6,0]) # when using the arm ends
    trimer1_13_rot1 = rotateFrame(trimer1_13_q1,init_coords)
    trimer1_24_q1 = VecsToQuat(trimer_v,[3*arm_end_dist/root6,3*arm_end_dist/root6,0])
    trimer1_24_rot1 = rotateFrame(trimer1_24_q1,init_coords)

    trimer2_13_q1 = VecsToQuat(trimer_v,[3*arm_end_dist/root6,-3*arm_end_dist/root6,0])
    trimer2_13_rot1 = rotateFrame(trimer2_13_q1,init_coords)
    trimer2_24_q1 = VecsToQuat(trimer_v,[-3*arm_end_dist/root6,3*arm_end_dist/root6,0])
    trimer2_24_rot1 = rotateFrame(trimer2_24_q1,init_coords)
    
    # 2nd roation
    trimer1_1_q2 = VecsToQuat([0,0,1],[1/root3,1/root3,1/root3])
    trimer1_1_rot2 = rotateFrame(trimer1_1_q2, trimer1_13_rot1)
    trimer1_3_q2 = VecsToQuat([0,0,1],[-1/root3,-1/root3,-1/root3])
    trimer1_3_rot2 = rotateFrame(trimer1_3_q2, trimer1_13_rot1)
    trimer1_2_q2 = VecsToQuat([0,0,1],[-1/root3,-1/root3,1/root3])
    trimer1_2_rot2 = rotateFrame(trimer1_2_q2, trimer1_24_rot1)
    trimer1_4_q2 = VecsToQuat([0,0,1],[1/root3,1/root3,-1/root3])
    trimer1_4_rot2 = rotateFrame(trimer1_4_q2, trimer1_24_rot1)

    trimer2_1_q2 = VecsToQuat([0,0,1],[-1/root3,1/root3,1/root3])
    trimer2_1_rot2 = rotateFrame(trimer2_1_q2, trimer2_13_rot1)
    trimer2_3_q2 = VecsToQuat([0,0,1],[1/root3,-1/root3,-1/root3])
    trimer2_3_rot2 = rotateFrame(trimer2_3_q2, trimer2_13_rot1)
    trimer2_2_q2 = VecsToQuat([0,0,1],[1/root3,-1/root3,1/root3])
    trimer2_2_rot2 = rotateFrame(trimer2_2_q2, trimer2_24_rot1)
    trimer2_4_q2 = VecsToQuat([0,0,1],[-1/root3,1/root3,-1/root3])
    trimer2_4_rot2 = rotateFrame(trimer2_4_q2, trimer2_24_rot1)
    
    # 8가지 trimer에 대한 translation 진행
    trimer1_1_coords = translation(trimer1_1_rot2,l/2,l/2,l/2)
    trimer1_3_coords = translation(trimer1_3_rot2,-l/2,-l/2,-l/2)
    trimer1_2_coords = translation(trimer1_2_rot2,-l/2,-l/2,l/2)
    trimer1_4_coords = translation(trimer1_4_rot2,l/2,l/2,-l/2)
    
    trimer2_1_coords = translation(trimer2_1_rot2,-l/2,l/2,l/2)  
    trimer2_3_coords = translation(trimer2_3_rot2,l/2,-l/2,-l/2)
    trimer2_2_coords = translation(trimer2_2_rot2,l/2,-l/2,l/2)
    trimer2_4_coords = translation(trimer2_4_rot2,-l/2,l/2,-l/2)

    save_new_coords(init_coords,"RK718_full_cube_origin.pdb")
    save_new_coords(trimer1_1_coords,"RK718_full_cube_1_1.pdb")
    save_new_coords(trimer1_2_coords,"RK718_full_cube_1_2.pdb")
    save_new_coords(trimer1_3_coords,"RK718_full_cube_1_3.pdb")
    save_new_coords(trimer1_4_coords,"RK718_full_cube_1_4.pdb")
    save_new_coords(trimer2_1_coords,"RK718_full_cube_2_1.pdb")
    save_new_coords(trimer2_2_coords,"RK718_full_cube_2_2.pdb")
    save_new_coords(trimer2_3_coords,"RK718_full_cube_2_3.pdb")
    save_new_coords(trimer2_4_coords,"RK718_full_cube_2_4.pdb")

cube_symm("C","307")


# In[ ]:




