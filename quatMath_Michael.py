import math
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

def rotation(position, rotate_quaternion):
    ref = [0,0,1]
    ## Setup rotation matrix
    q0,q1,q2,q3 = rotate_quaternion[0],rotate_quaternion[1],rotate_quaternion[2],rotate_quaternion[3] 
    r00 = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
    r01 = 2.0 * (q1 * q2 - q0 * q3);
    r02 = 2.0 * (q0 * q2 + q1 * q3);
    r10 = 2.0 * (q0 * q3 + q1 * q2);
    r11 = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
    r12 = 2.0 * (q2 * q3 - q0 * q1);
    r20 = 2.0 * (q1 * q3 - q0 * q2);
    r21 = 2.0 * (q0 * q1 + q2 * q3);
    r22 = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
    
    ## Rotate reference vector
    x = r00 * position[0] + r01 * position[1] + r02 * position[2]
    y = r10 * position[0] + r11 * position[1] + r12 * position[2]
    z = r20 * position[0] + r21 * position[1] + r22 * position[2]
    
    return np.array([x,y,z])

def calAngle(v1,v2):
    # Angle between v1 and v2
    cosine_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.arccos(cosine_angle) # Radian
    angle = np.degrees(angle) # Degree
    return angle

def calDist(v1,v2):
    return math.sqrt((v2[0]-v1[0])**2 + (v2[1]-v1[1])**2 + (v2[2]-v1[2])**2)

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

