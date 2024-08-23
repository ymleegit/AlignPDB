from __future__ import division
from __future__ import print_function
import numpy as np
import math
import random
import numpy.testing as npt

## Returns the angle required to rotate the orientation prescribed by quatA into the orientation prescribed by quatB
## i.e., if quatB = R*quatA, then R = quatB*quatA^dag. Angle associated with R is contained in its 0th component.
# input quats are arrays of four numbers: [q0, q1, q2, q3]
def computeSeparationAngle(quatA, quatB):
    R = quatMultiply(quatB, quatConjugate(quatA))
    # this should be impossible but due to round-off errors we must be careful
    # default atol (absolute tolerance) value for np.isclose is 1e-8
    if np.isclose(R[0], 1.0, atol=1e-7, rtol=0.0):
        R[0] = 1.0
    if np.isclose(R[0], -1.0, atol=1e-7, rtol=0.0):
        R[0] = -1.0
    # the range of acos is [0,pi], so we are fine.
    theta = 2.*math.acos(R[0])

    return theta

## Returns the quaternion that is EQUIVALENT to quat that has minimal angle with ref_quat.
## Equivalent quaternions are determined according to the list of equiv_quats provided by
## the user. They must encode all rotational symmetry of the particle whose orientation quat prescribes.
# input quats are arrays of four numbers: [q0, q1, q2, q3]
def getMinEquivOrientation(ref_quat, quat, equiv_quats):
    # start with the quaternion before it's been rotated by any equivalent quaternions.
    min_quat = quat
    min_angle = computeSeparationAngle(quat, ref_quat)

    # loop through all equivalent quaternions and see if they have smaller angles with ref_quat.
    for q in equiv_quats:

        test_quat = quatMultiply(quat, q)
        # we must include this, since q and -q effect the same rotation.
        mirror_test_quat = -1.*test_quat

        for eq in [test_quat, mirror_test_quat]:
            angle = computeSeparationAngle(eq, ref_quat)

            if angle < min_angle:
                min_angle = angle
                min_quat = eq

    return min_quat, min_angle

## Returns a random quaternion culled from a uniform distribution on the surface of a 3-sphere.
## Uses the MARSAGLIA (1972) method (a la hoomd)
## NOTE THAT generating a random rotation via a random angle about a random axis of rotation is INCORRECT.
## See K. Shoemake, "Uniform Random Rotations," 1992, for a nice explanation for this.
# output quat is an array of four numbers: [q0, q1, q2, q3]
def quatRandom():
    # random.uniform(a,b) gives number in [a,b]
    v1 = random.uniform(-1,1)
    v2 = random.uniform(-1,1)
    v3 = random.uniform(-1,1)
    v4 = random.uniform(-1,1)

    s1 = v1*v1 + v2*v2
    s2 = v3*v3 + v4*v4

    while (s1 >= 1.):
        v1 = random.uniform(-1,1)
        v2 = random.uniform(-1,1)
        s1 = v1*v1 + v2*v2

    while (s2 >= 1. or s2 == 0.):
        v3 = random.uniform(-1,1)
        v4 = random.uniform(-1,1)
        s2 = v3*v3 + v4*v4

    s3 = np.sqrt((1.-s1)/s2)

    return np.array([v1, v2, v3*s3, v4*s3])

## Returns the conjugate of an input quaternion.
# input quats are arrays of four numbers: [q0, q1, q2, q3]
def quatConjugate(quat):
    # quaternion scalar
    qs = quat[0]
    # quaternion vector
    qv = np.array([quat[1], quat[2], quat[3]])
    # invert the vector
    qv  = -1.0*qv

    return np.array([qs, qv[0], qv[1], qv[2]])

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

## Returns the rotation of a vector by a given quaternion.
## Shamelessly from Greg's quatvec function in functors.h in his MC code
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
    vecarr=np.asarray(veclist)
    new_veclist=[]

    for vec in vecarr:
        newvec = quatVec(quat, vec)
        new_veclist.append(newvec)

    new_veclist=np.asarray(new_veclist)

    return new_veclist

## Returns the scalar invariant combination S of quaternions quat1 and quat2. Takes into account particle symmetry, if that is provided by the user,
## by calculating S for all possible pairs of equivalent quaternions, and returning that list of S values.
# input quats are arrays of four numbers: [q0, q1, q2, q3]
# equiv_quats1: list of quaternions that encode the rotational symmetry of particle 1. Each equiv_quat1 represents
# a rotation that leaves particle 1 unchanged.
# equiv_quats2: list of quaternions that encode the rotational symmetry of particle 2. Each equiv_quat2 represents
# a rotation that leaves particle 2 unchanged.
def computeSymmS(quat1, quat2, equiv_quats1=[np.array([1,0,0,0])], equiv_quats2=[np.array([1,0,0,0])],speed=True):
    S_arr = []
    if speed == True:
        # loop through all equivalent quaternions for quat1
        for eq1 in equiv_quats1:
            tq1 = quatMultiply(quat1, eq1)
            # loop through all equivalent quaternions for quat2
            for eq2 in equiv_quats2:
                tq2  = quatMultiply(quat2, eq2)
                # scalar invariant S
                S = computeS(tq1, tq2)
                # eq1 and -eq1 effect the same rotation, and eq2 and -eq2 effect the same rotation.
                # thus, -tq1 and -tq2 are equally valid equivalent quaternions.
                # S(p,-q) = S(-p,q) = -S(p,q), so we must include -S in our list. [note that S(-p,-q) = S(p,q)].
                Sneg = -1.*S
                S_arr.append(S)
                S_arr.append(Sneg)
    else:
        # loop through all equivalent quaternions for quat1.
        for eq1 in equiv_quats1:
            tq1 = quatMultiply(quat1, eq1)
            # we must include this, since q and -q effect the same rotation.
            mtq1 = -1.*tq1
            for q1 in [tq1, mtq1]:
                for eq2 in equiv_quats2:
                    tq2  = quatMultiply(quat2, eq2)
                    mtq2 = -1.*tq2
                    for q2 in [tq2, mtq2]:
                        # scalar invariant S
                        S = computeS(q1, q2)
                        S_arr.append(S)

    return S_arr

## Returns the scalar invariant combination S of quaternions quat1 and quat2, as defined in G. van Anders et al., PNAS 111, E4812 (2014), supplementary info.
# input quats are arrays of four numbers: [q0, q1, q2, q3]
def computeS(q1, q2):
    # scalar invariant S
    S = 0.5*(quatMultiply(q1, quatConjugate(q2)) + quatMultiply(q2, quatConjugate(q1)))

    # assert that the vector part of S is zero
    npt.assert_almost_equal(np.array([S[1], S[2], S[3]]), np.zeros(3), decimal=6, err_msg="S is not a purely scalar quantity!")

    return S[0]

## Returns the scalar invariant combination U of quaternions quat1 and quat2. Takes into account particle symmetry, if that is provided by the user,
## by calculating U for all possible pairs of equivalent quaternions, and returning that list of U values.
# input quats are arrays of four numbers: [q0, q1, q2, q3]
# equiv_quats1: list of quaternions that encode the rotational symmetry of particle 1. Each equiv_quat1 represents
# a rotation that leaves particle 1 unchanged.
# equiv_quats2: list of quaternions that encode the rotational symmetry of particle 2. Each equiv_quat2 represents
# a rotation that leaves particle 2 unchanged.
def computeSymmU(quat1, quat2, equiv_quats1=[np.array([1,0,0,0])], equiv_quats2=[np.array([1,0,0,0])]):
    U_arr = []
    # loop through all equivalent quaternions for quat1.
    for eq1 in equiv_quats1:
        tq1 = quatMultiply(quat1, eq1)
        # we must include this, since q and -q effect the same rotation.
        mtq1 = -1.*tq1
        for q1 in [tq1, mtq1]:
            for eq2 in equiv_quats2:
                tq2  = quatMultiply(quat2, eq2)
                mtq2 = -1.*tq2
                for q2 in [tq2, mtq2]:
                    # scalar invariant S
                    U = computeU(q1, q2)
                    U_arr.append(U)

    return U_arr

## Returns the scalar invariant combination U of quaternions quat1 and quat2, as defined in G. van Anders et al., PNAS 111, E4812 (2014), supplementary info.
# input quats are arrays of four numbers: [q0, q1, q2, q3]
def computeU(q1, q2):
    # scalar invariant U
    U = 0.25*(quatMultiply(q1, q2) + quatMultiply(q2, q1) + quatMultiply(quatConjugate(q1), quatConjugate(q2)) + quatMultiply(quatConjugate(q2), quatConjugate(q1)))

    # assert that the vector part of U is zero
    npt.assert_almost_equal(np.array([U[1], U[2], U[3]]), np.zeros(3), decimal=6, err_msg="U is not a purely scalar quantity!")

    return U[0]

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

## gets the quaternion that rotates vec set v1 to vec set v2
# v1 and v2 MUST have shape (3,N)
# v1 and v2 MUST be wrt the origin.
def VecSetsToQuat(v1, v2):
    # get the rotation matrix
    R = KabschAlgorithm(v1, v2)
    # check that this matrix results in essentially ZERO RMSD between rotated vector sets
    v1_rot = np.dot(R, v1)
    for vec1, vec2 in zip(np.transpose(v1_rot), np.transpose(v2)):
        npt.assert_almost_equal(vec1, vec2, decimal=5)
    # turn the rotation matrix into a quaternion
    q = RotMatToQuat(R)

    return q

## turns a rotation matrix into a quaternion in a stable manner
# see: Animating Rotation with Quaternion Curves, Ken Shoemake, 1985.
# R MUST be a (3x3) numpy array that is orthogonal etc.
def RotMatToQuat(R):
    wsq = 0.25*(1+R[0,0]+R[1,1]+R[2,2])
    if wsq > np.finfo(np.single).eps:
        w = np.sqrt(wsq)
        x = (R[1,2]-R[2,1])/(4*w)
        y = (R[2,0]-R[0,2])/(4*w)
        z = (R[0,1]-R[1,0])/(4*w)
    else:
        w = 0.
        xsq = -0.5*(R[1,1]+R[2,2])
        if xsq > np.finfo(np.single).eps:
            x = np.sqrt(xsq)
            y = R[0,1]/(2*x)
            z = R[0,2]/(2*x)
        else:
            x = 0.
            ysq = 0.5*(1.-R[2,2])
            if ysq > np.finfo(np.single).eps:
                y = np.sqrt(ysq)
                z = R[1,2]/(2*y)
            else:
                y = 0.
                z = 1.

    # the above actually finds the adjoint (inverse) of the quaternion, the way I've defined it.
    # so return the inverse of this (the inverse of the inverse)
    return np.array([w,-x,-y,-z])

## gets the rotation matrix that rotates v1 to v2 and results in minimal RMSD between the vector sets
# v1 and v2 MUST have shape (3,N)
# v1 and v2 MUST be wrt the origin.
# https://en.wikipedia.org/wiki/Kabsch_algorithm. See also Paul's script brute_force.h, in the registration freud branch.
def KabschAlgorithm(v1, v2):
    assert(len(v1) == len(v2) == 3)
    # calculate the covariance matrix
    A = np.dot(v1, np.transpose(v2))
    # calculate the singular value decomposition of A
    V, s, W = np.linalg.svd(A)
    # this is dumb. linalg.svd actually returns the transpose of W, rather than W itself.
    W = np.transpose(W)
    # ensure right-handed coordinate system
    d = np.sign(np.linalg.det(W)*np.linalg.det(V))
    D = np.array([[1,0,0],[0,1,0],[0,0,d]])
    # get the rotation matrix
    R = np.dot(D, np.transpose(V))
    R = np.dot(W, R)

    return R
