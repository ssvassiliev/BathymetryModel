#!/usr/bin/python

import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import cm
import csv, math, numpy, sys, os

# Required input:
WorkDir = os.getcwd()+"/"
verticesFile = WorkDir+'Opinicon/Output/vertices.csv'
facesFile = WorkDir+'Opinicon/Output/faces.csv'

# Options
#--- Scale depth
zsc = -1.0

print "<<< Reading vertices >>>"
# read vertices
x=[]; y=[]; z=[]
try:
  with open(verticesFile) as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
      x.append(float(row[0]))
      y.append(float(row[1]))
      z.append(float(row[2]))
except IOError:
    print 'Error: file',file,'not found'
    raise SystemExit
n1=len(x)
print "... Number of vertices:", n1

for i in range(n1):
    z[i] *= zsc

print "<<< Reading triangles >>>"
# read triangles
triangles=[]
try:
    with open(facesFile) as csvDataFile:
        csvReader = csv.reader(csvDataFile)
        for row in csvReader:
            triangles.append([int(row[0]),int(row[1]),int(row[2])])
except IOError:
      print 'Error: file',file,'not found'
      raise SystemExit
n1=len(triangles)
print "... Number of triangles:", n1

plt.figure()
plt.rcParams['axes.facecolor'] = 'brown'
plt.gca().set_aspect('equal')
#plt.tricontour(x, y, triangles, z, 24, cmap=cm.ocean, linewidths=1)
plt.tricontourf(x, y, triangles, z, 80, cmap=cm.ocean)
c1=plt.colorbar()
c1.set_label('Depth, m')
plt.grid(c='y', ls='-', alpha=0.3)
plt.minorticks_on()
plt.title('L.Opinicon bathymetry')
plt.xlabel('Easting')
plt.ylabel('Northing')
# Set x,y limits or comment out to draw the whole map 
#plt.xlim(394800,395800) 
#plt.ylim(4934500,4935500)
#plt.xlim(391320, 395935)
#plt.ylim(4931861, 4936171)
plt.tight_layout()
plt.savefig('opinicon_bathymetry.png', dpi = 600)

plt.figure()
plt.rcParams['axes.facecolor'] = 'white'
plt.gca().set_aspect('equal')
plt.triplot(x, y, triangles, lw=0.5, color='black')
# Interactive plot
plt.show()






