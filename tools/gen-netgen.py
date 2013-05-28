#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Generate netgen mesh from python function """

# make 2/1 = 0.5, instead of 0
from __future__ import division

import numpy as np
import matplotlib.pyplot as plt


filename = 'test' # .in2d will be added
##### Variables

n=100	# Number of points in 1 sinus
L=1.5	# Wavelenght of the roughtness
T=10.0	# Number of periods 
a=0.15	# Hight of the Sinus
s=15	# Heaviside coefficient
GradingFactor=2

## Geometrical parameters
H0=1.5	# Hight of the WG
Lsin=L*T	# Lenght of the roughtness
Lb=1.5	# Lenght before the roughtness
La=1.5	# Lenght after the roughtness
l=1.5	# Thickness of PML1
H1=3.0	# Hight of the substrat
H2=6.0	# Hight of the PML1

# Domains :
dwg=1	# WG
dst=2	# substrat at the top
dsb=3	# substrat at the bottom
dpml1t=4	# PML1 at the top
dpml1b=5	# PML1 at the bottom
dpml2l=6	# PML2 left
dpml2r=7	# PML2 right
dout=0	# outside

# Boundary conditions
bcwgs=1
bcwgpml2=2
bcspml1=3
bcpml12=4
bcout=5

## Shape of the roughtness
x=np.linspace(-0.1, Lsin+0.1, n)		
h=1/(1+np.exp(-s*(x))) - 1/(1+np.exp(-s*(x-15)))	# Heaviside Function
y=a*np.sin(x*2*np.pi/L) * h			# Function for the roughtness

##### Functions

def Point(P,n,x,y):
    """ Add point to array. Create new array if array is empty """
    if np.any(P):
        return np.vstack((P,[n,x,y]))
    else:
        return np.array([n,x,y])


## Rectangles
# (X0,Y0) are the coordinates for the middle of the left side of the rectangle
# L : lenght of the rectangle
# H : hight of the rectangle

def RectPoints(P,n,R,X0,Y0,L,H):	# R=1 for the first rectangle, R=2 for the second etc.
    """ 4 points defined for a rectangle """
    P=Point(P,(2*n+4*R+1),(X0+L),(Y0-H/2))
    P=Point(P,(2*n+4*R+2),(X0+L),(Y0+H/2))
    P=Point(P,(2*n+4*R+3),X0,(Y0+H/2))
    P=Point(P,(2*n+4*R+4),X0,(Y0-H/2))
    return P

def WGcornersPoints(P,n,Lsin,La,Lb,H0):
    """ 4 points defined for the corners of the WG """
    P=Point(P,(2*n+1),(Lsin+La),(-H0/2.))
    P=Point(P,(2*n+2),(Lsin+La),(+H0/2.))
    P=Point(P,(2*n+3),(-Lb),(H0/2.))
    P=Point(P,(2*n+4),(-Lb),(-H0/2.))
    return P

#def Segment(dl,dr,p1,p2,bc):
#    """ 1 segment defined """
#    f.write('%i \t %i \t 2 \t %.i \t %.i \t -bc=%d \n' % (dl,dr,p1,p2,bc))
    
## Segments : A is an array, each line defines a segment : (dl, dr, Nbr of pts = 2, P1, P2, bc)

# The function "add" adds a line in A
def add(A,dl,dr,p1,p2,bc):
    if np.any(A):
        return np.vstack((A,[dl,dr,2,p1,p2,bc]))
    else:
        return np.array([dl,dr,2,p1,p2,bc])

# Polygon is used to draw the limits of the PML1 and of the substrat.
def Polygon(A,p1,p2,p3,p4,p5,p6,p7,p8, din1,din2,din3,din4,dout1,dout2,dout3,dout4, bc1,bc2,bc3,bc4):
    """ Polygon defined, with 8 segments """
    A=add(A,din1,dout1,p1,p2,bc1)
    A=add(A,din2,dout2,p2,p3,bc2)
    A=add(A,din3,dout3,p3,p4,bc3)
    A=add(A,din3,dout3,p4,p5,bc3)
    A=add(A,din3,dout3,p5,p6,bc3)
    A=add(A,din4,dout4,p6,p7,bc4)
    A=add(A,din1,dout1,p7,p8,bc1)
    A=add(A,din1,dout1,p8,p1,bc1)
    return A

# Rectangle is used to draw the 2 PML2, at the left and at the right of the WG
def Rectangle(A,p1,p2,p3,p4, din,d1,d2,d3,d4, bc1,bc2,bc3,bc4):	
    """ 4 segments defined to make a rectangle """
    A=add(A,din,d1,p1,p2,bc1)
    A=add(A,din,d2,p2,p3,bc2)
    A=add(A,din,d3,p3,p4,bc3)
    A=add(A,din,d4,p4,p1,bc4)
    return A

# To draw the corners of the WG
def WGcorners(A,n,dinwg,doutrwg,douttwg,doutlwg,doutbwg, bcwgs,bcwg2):
    """ 6 segments defined to close the WG """
    A=add(A,dinwg,doutbwg,(2*n),(2*n+1),bcwgs)
    A=add(A,dinwg,doutrwg,(2*n+1),(2*n+2),bcwg2)
    A=add(A,dinwg,douttwg,(2*n+2),n,bcwgs)
    A=add(A,dinwg,douttwg,1,(2*n+3),bcwgs)
    A=add(A,dinwg,doutlwg,(2*n+3),(2*n+4),bcwg2)
    A=add(A,dinwg,doutbwg,(2*n+4),(n+1),bcwgs)
    return A


############## Generate Points
P = []

## sinus at the top	
for i in range(len(x)):
    P = Point(P,i+1,x[i],y[i]+H0/2)

## sinus at the bottom	
for i in range(len(x)):
    P = Point(P,len(x)+i+1,x[i],y[i]-H0/2)

## Points for the corners of the WG 
P = WGcornersPoints(P,n,Lsin,La,Lb,H0)

## Points for the rectangles : here are defined all the points needed to draw the PML1, PML2, and the Substrat.
P = RectPoints(P,n,1,(-Lb),0,(Lb+Lsin+La),H1)	
P = RectPoints(P,n,2,(-Lb-l),0,(Lb+Lsin+La+2*l),H0)	
P = RectPoints(P,n,3,(-Lb-l),0,(Lb+Lsin+La+2*l),H2)

## Probe			(Useless if we don't want to draw a 'probe')
#Re=3	# Number of rectangles
#Xp=16.0	# Position of the probe
#Hp=1.5	# Hight of the probe
#Point((2*n+4*Re+5),Xp,(Hp/2))
#Point((2*n+4*Re+6),Xp,(-Hp/2))


######### Generate Segments
# Definition of the array, which will contain all the segments
# A=["#dl","dr",2,"p1","p2","-bc"]
A = []

## Sin at the top
indh=np.arange(1,len(x))
for ih in indh:
    A=add(A,dst,dwg,ih,(ih+1),bcwgs)

## Sin at the bottom
indh=np.arange(len(x)+1,2*len(x))
for ih in indh:
    A=add(A,dwg,dsb,ih,(ih+1),bcwgs)

## Corners of the WG
A=WGcorners(A,n,dwg,dpml2r,dst,dpml2l,dsb, bcwgs,bcwgpml2)

## First rectangle : PML2 left
A=Rectangle(A,(2*n+4),(2*n+3),(2*n+11),(2*n+12), dpml2l,dwg,dpml1t,dout,dpml1b, bcwgpml2,bcpml12,bcout,bcpml12)

## Second rectangle : PML2 right
A=Rectangle(A,(2*n+9),(2*n+10),(2*n+2),(2*n+1),dpml2r,dout,dpml1t,dwg,dpml1b, bcout,bcpml12,bcwgpml2,bcpml12)

# Substrat
#Rectangle((2*n+5),(2*n+6),(2*n+7),(2*n+8),dst,dpml1t,dpml1t,dpml1b,dpml1b, bcout,bcpml12,bcwgpml2,bcpml12)
A=Polygon(A,(2*n+5),(2*n+1),(2*n+2),(2*n+6),(2*n+7),(2*n+3),(2*n+4),(2*n+8), dsb,dwg,dst,dwg, dpml1b,dpml2r,dpml1t,dpml2l, bcspml1,bcwgpml2,bcspml1,bcwgpml2)	

# PML
#Rectangle((2*n+13),(2*n+14),(2*n+15),(2*n+16),dpml1t,dout,dout,dout,dout, bcout,bcpml12,bcwgpml2,bcpml12)
A=Polygon(A,(2*n+13),(2*n+9),(2*n+10),(2*n+14),(2*n+15),(2*n+11),(2*n+12),(2*n+16), dpml1b,dpml2r,dpml1t,dpml2l, dout,dout,dout,dout, bcout,bcout,bcout,bcout)

## Probe	 		(Useless if we don't want to draw a 'probe')
#Segment(dinwg,dinwg,(2*n+4*Re+5),(2*n+4*Re+6),bcwg)

## Loop : we have to delete all double segments, so as to each line is defined only once.

deleteme = np.array([])

for i in range(len(A)):
    if i not in deleteme:
        var1=np.where(
                ( (A[:,3]==A[i,3]) & (A[:,4]==A[i,4]) ) 
                | 
                ( (A[:,3]==A[i,4]) & (A[:,4]==A[i,3]) )
                )[0]
        deleteme = np.append(deleteme,var1[1:])

A=np.delete(A,deleteme,0)


####### We write everything into the file

f=open(filename + '.in2d','w')	# open the text file
f.write('splinecurves2dv2 \n%.i \npoints \n' % GradingFactor )	# these first lines are nessecary

####### Points
#(Pt Nbr, x, y)

for i in range(len(P)):
    f.write("%i \t %f \t %f \n" % (int(P[i,0]),float(P[i,1]),float(P[i,2])))

###### Segments
#(domain left    domain right    Nbr of Pts    point1    point2    bc)

f.write('\nsegments \n')
for i in range(len(A)):
    f.write("%i \t %i \t %.i \t %.i \t %.i \t -bc=%.i \n" % (int(A[i,0]),int(A[i,1]),int(A[i,2]),int(A[i,3]),int(A[i,4]),int(A[i,5])))

#### Materials
f.write('\nmaterials \n')
f.write('1 \t domain1 \t -maxh=1 \n')
f.write('2 \t domain2 \t -maxh=2 \n')
f.write('3 \t domain3 \t -maxh=2 \n')
f.write('4 \t domain4 \t -maxh=2 \n')
f.write('5 \t domain5 \t -maxh=2 \n')
f.write('6 \t domain6 \t -maxh=2 \n')
f.write('7 \t domain7 \t -maxh=2 \n')


f.close()    
