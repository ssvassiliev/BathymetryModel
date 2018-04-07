#!/usr/bin/python

from stl import mesh
import csv, math, numpy, sys, os
import triangle
import triangle.plot
print "\n************** Triangulate 2 ***************"
#---------------------------------------------------------
# This round of triangulation is designed to extend lake
# bounds to zero lines, extrude zero lines to enhance printed
# walls, center and scale coordinates 
#---------------------------------------------------------
# Make holes in place of the islands?
makeHoles = True
extrudeWalls = False
# Size of the largest mesh dimension
size = 20 
WorkDir = os.getcwd()+"/"
verticesFile = WorkDir+'Opinicon/Output/vertices.csv'
vertices_extFile = WorkDir+'Opinicon/Output/vertices_ext.csv'
holesFile = WorkDir+'Opinicon/Output/holes.csv'
segmentsFile = WorkDir+'Opinicon/Output/segments.csv'

print "------ Variables -------"
print "makeHoles =", makeHoles
print "extrudeWalls =", True
print "output size =", size
print "------------------------"

print "<<< Reading holes >>>"
# read holes
holes=[]
try:
  with open(holesFile) as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
      holes.append([float(row[0]),float(row[1])])
except IOError:
  print 'Error: file',holesFile,'not found'
  raise SystemExit   
ns=len(holes)
HOLES = numpy.ndarray(shape = (ns,2), dtype = float)
for i in range(ns):
   HOLES[i,0] = holes[i][0]; HOLES[i,1] = holes[i][1]; 

print "<<< Reading vertices >>>"
# read vertices
vertices=[]
try:
  with open(verticesFile) as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
      vertices.append([float(row[0]),float(row[1]),float(row[2])])
except IOError:
  print 'Error: file',file,'not found'
  raise SystemExit
n1=len(vertices)
print "... Number of vertices:", n1

if extrudeWalls:
   # read extruded wall (for printing)
   print "<<< Reading extuded vertices >>>"
   try:
     with open(vertices_extFile) as csvDataFile:
       csvReader = csv.reader(csvDataFile)
       for row in csvReader:
         vertices.append([float(row[0]),float(row[1]),float(row[2])])
   except IOError:
     print 'Error: file',file,'not found'
     raise SystemExit
ns=len(vertices)
print "... Number of vertices:", ns

XY = numpy.ndarray(shape = (ns,2), dtype = float)
for i in range(ns):
   XY[i,0] = vertices[i][0]; XY[i,1] = vertices[i][1];

# read segments
print "<<< Reading segments >>>"
segments=[]
try:
  with open(segmentsFile) as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
      segments.append([int(row[0]),int(row[1])])
except IOError:
  print 'Error: file',file,'not found'
  raise SystemExit

# extrude wall vertices
print "<<< Extruding outer wall for printing >>>"
for i in range(n1,ns-1):
   segments.append([i,i+1])
segments.append([i+1,n1])
ns=len(segments)
SEGM = numpy.ndarray(shape = (ns,2), dtype = int)
for i in range(ns):
   SEGM[i,0] = segments[i][0]; SEGM[i,1] = segments[i][1]; 

print "<<< Generating meshes >>>"   
#----------------------------------------------------
# triangulate
if makeHoles:
   A = dict(vertices=XY, segments=SEGM, holes=HOLES)
else:   
   A = dict(vertices=XY, segments=SEGM)
B = triangle.triangulate(A,'pq10')

# prepare 3D verices by adding z coordinate
na=len(B['vertices'])
vrtb = numpy.ndarray(shape = (na,3), dtype = float)
vrtt = numpy.ndarray(shape = (na,3), dtype = float)

#old vertices
for i in range(n1):
   vrtb[i,0] = vrtt[i,0] = B['vertices'][i][0]
   vrtb[i,1] = vrtt[i,1] = B['vertices'][i][1]
   vrtb[i,2] = vertices[i][2]; vrtt[i,2] = 0.0
#new vertices
for i in range(n1,na):
   vrtb[i,0] = vrtt[i,0] = B['vertices'][i][0]
   vrtb[i,1] = vrtt[i,1] = B['vertices'][i][1]
   vrtb[i,2] = vrtt[i,2] = 0.0
   
# <<<<<<<<<<<<<    Center coordinates    >>>>>>>>>>>>>>
#----------------------------------------------------------
print "<<< Centering meshes >>>"
center = numpy.mean(vrtb, 0)
center[2]=0.0
print "... Coordinates of the geometric center:\n","...", center
vrtb = vrtb - center
vrtt = vrtt - center
HOLES = HOLES - center[[0,1]]
# <<<<<<<<<<<<<    Scale coordinates    >>>>>>>>>>>>>>
#----------------------------------------------------------
xSize=max(vrtb[:,0])-min(vrtb[:,0])
ySize=max(vrtb[:,1])-min(vrtb[:,1])
print "... Mesh size:",xSize,"x",ySize
print "<<< Scaling meshes >>>"
vrtb *= size/max(xSize,ySize)
vrtt *= size/max(xSize,ySize)

# the faces (triangles)
faces = B['triangles']
# Create meshes
top_msh = mesh.Mesh(numpy.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
bottom_msh = mesh.Mesh(numpy.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        bottom_msh.vectors[i][j] = vrtb[f[j],:]
        top_msh.vectors[i][j] = vrtt[f[j],:]
# Write meshes to files
bottom_msh.save('bottom_mesh.stl')
top_msh.save('top_mesh.stl')
