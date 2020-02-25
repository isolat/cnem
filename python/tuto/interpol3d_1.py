##
## Copyright (C) 2003-2011
## Lounes ILLOUL (illoul_lounes@yahoo.fr)
## Philippe LORONG (philippe.lorong@ensam.eu)
## Arts et Metiers ParisTech, Paris, France
##
## $Revision: 25 $ - $Date: 2012-02-22 17:36:30 +0100 (mer., 22 fevr. 2012) $

import shutil
import re
import math
import cPickle
import numpy as np
import scipy as sp

import sys
sys.path.append('../bin')

import naturalInterpolant as NT

import time

t_ini = time.time()

#-------------------------------------------------------------------------------
# read nodes and boundary :
print '\n--------------------------------------------------------------------------------\n'
print 'creat dummy point for test...\n'

X,Y,Z=np.mgrid[-0.5:0.5:3j,-0.5:0.5:3j,-0.5:0.5:3j]
XYZ_Noeud=np.zeros((X.reshape(-1).shape[0],3),dtype=np.float64,order='C')
XYZ_Noeud[:,0],XYZ_Noeud[:,1],XYZ_Noeud[:,2]=X.reshape(-1),Y.reshape(-1),Z.reshape(-1)
#print 'nb node : ',X.reshape(-1).shape[0]
XYZ_Noeud=XYZ_Noeud.reshape(-1)

#-------------------------------------------------------------------------------
# creat dummy point location for interpolation evaluation

#XYZ_Point=list((np.random.random(1000000*3)-0.5)*0.15)

X,Y,Z=np.mgrid[-0.45:0.45:10j,-0.45:0.45:10j,-0.45:0.45:10j]
XYZ_Point=np.zeros((X.reshape(-1).shape[0],3),dtype=np.float64,order='C')
XYZ_Point[:,0],XYZ_Point[:,1],XYZ_Point[:,2]=X.reshape(-1),Y.reshape(-1),Z.reshape(-1)
#print 'nb point : ',X.reshape(-1).shape[0]
XYZ_Point=XYZ_Point.reshape(-1)

#-------------------------------------------------------------------------------
# creat the interpolator
print '\n--------------------------------------------------------------------------------\n'
print 'creat the interpolator...\n'

Type_FF=0
Fxyz=NT.interpol(XYZ_Noeud,[],XYZ_Point,Type_FF)

## if we don't enter the boundary tri IN_Tri_Ini, the convex hull of the nodes
##    will be taken as boundary:
#interpol=NT.interpol(XYZ_Noeud,[],XYZ_Point,Type_FF)

#----------------------------------------------------------------------------
# interpolate variable filds : Var
print '\n--------------------------------------------------------------------------------\n'
print 'interpolate variable filds...\n'

# test : Var = XYZ_Noeud ==> interpolated Var on XYZ_Point = XYZ_Point

Var=np.array(XYZ_Noeud)
Var=Var.reshape((Var.shape[0]/3,3))

#interpolate Var :
#Var_Int=Fxyz.interpolate(Var)

Mat_Int=Fxyz.mat_interpol_glob()
Var_Int=Mat_Int*Var

elapsed = time.time() - t_ini
print 'elapsed : ',elapsed

XYZ_Point=np.array(XYZ_Point)
XYZ_Point=XYZ_Point.reshape((XYZ_Point.shape[0]/3,3))

Dif=Var_Int-XYZ_Point

maxdif=0.

for i in xrange(len(Dif)):
    if Fxyz.In_Out[i] : # point interpolation in the domain 
        maxdif_i=max(abs(Dif[i]))
        if maxdif_i>maxdif:
            maxdif=maxdif_i
        
print 'err max : ', maxdif

#-------------------------------------------------------------------------------

#py_cnem3d.interpol.in_use()
