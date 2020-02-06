#!/usr/bin/env python3

'''Unroll QI = R G Rt, where QI and G are both traceless quadrupoles, and R
rotates a vector from global frame to quasi-internal frame.'''

from sympy import symbols, Matrix
from sympy.codegen.rewriting import create_expand_pow_optimization

def get_qi():
    '''QI = R G Rt'''

    g00,g11,g22 = symbols('a1xx,a2yy,a3zz')
    g01,g02,g12 = symbols('a4xy,a5xz,a6yz')
    r00,r01,r02 = symbols('r00,r01,r02')
    r10,r11,r12 = symbols('r10,r11,r12')
    r20,r21,r22 = symbols('r20,r21,r22')
    gmat = Matrix([[g00,g01,g02],[g01,g11,g12],[g02,g12,g22]])
    rot = Matrix([[r00,r01,r02],[r10,r11,r12],[r20,r21,r22]])
    qm = Matrix([[0,0,0],[0,0,0],[0,0,0]])
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for m in range(3):
                    qm[i,j] += rot[i,k]*gmat[k,m]*rot[j,m]
    return qm

if __name__ == '__main__':
    expand_pow = create_expand_pow_optimization(2)
    qi = get_qi()
    print('qixx = {};'.format(expand_pow(qi[0,0])))
    print('qiyy = {};'.format(expand_pow(qi[1,1])))
    print('// qizz = {};'.format(expand_pow(qi[2,2])))
    print('qixy = {};'.format(expand_pow(qi[0,1])))
    print('qixz = {};'.format(expand_pow(qi[0,2])))
    print('qiyz = {};'.format(expand_pow(qi[1,2])))
    print('qizz = -(qixx + qiyy);')
