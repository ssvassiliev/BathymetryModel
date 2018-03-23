#!/usr/bin/python
import scipy, numpy
from scipy.spatial.distance import cdist  
from mergepoints import merge_points, sqdistance
from stl import mesh
import shapefile, fiona, csv, math, numpy, scipy, os, sys

#---------------------------------------------------------
# <<<<<< Triangulate Planar Straight Line Graph >>>>>>>>>
#---------------------------------------------------------
import triangle
import triangle.plot

print "<<< Constrained conforming Delaunay triangulation >>>\n"   
#----------------------------------------------------
# <<<<<<< read vertices, segments, holes >>>>>>>

# read segments
segments=[]
file = 'segments.csv'
try:
  with open(file) as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
      segments.append([int(row[0]),int(row[1])])
except IOError:
  pass
ns=len(segments)
SEGM = numpy.ndarray(shape = (ns,2), dtype = int)
for i in range(ns):
   SEGM[i,0] = segments[i][0]; SEGM[i,1] = segments[i][1]; 

# read holes
holes=[]
file = 'holes.csv'
try:
  with open(file) as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
      holes.append([float(row[0]),float(row[1])])
except IOError:
  pass   
ns=len(holes)
HOLES = numpy.ndarray(shape = (ns,2), dtype = float)
for i in range(ns):
   HOLES[i,0] = holes[i][0]; HOLES[i,1] = holes[i][1]; 


# read vertices
vertices=[]
file = 'vertices.csv'
try:
  with open(file) as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
      vertices.append([float(row[0]),float(row[1]),float(row[2])])
except IOError:
  pass   
ns=len(vertices)  
XY = numpy.ndarray(shape = (ns,2), dtype = float)
for i in range(ns):
   XY[i,0] = vertices[i][0]; XY[i,1] = vertices[i][1]; 
   
A = dict(vertices=XY, segments=SEGM, holes=HOLES)
B = triangle.triangulate(A,'pq0')

ns=len(B['vertices'])
vrtb = numpy.ndarray(shape = (ns,3), dtype = float)
vrtt = numpy.ndarray(shape = (ns,3), dtype = float)
for i in range(ns-1):
   vrtb[i,0] = vrtt[i,0] = B['vertices'][i][0]
   vrtb[i,1] = vrtt[i,1] = B['vertices'][i][1]
   vrtb[i,2] = vertices[i][2]; vrtt[i,2] = 0.0

   vrtb[ns-1,0] = vrtt[ns-1,0] = B['vertices'][ns-1][0]
   vrtb[ns-1,1] = vrtt[ns-1,1] = B['vertices'][ns-1][1]
   vrtb[ns-1,2] = 0.0; vrtt[ns-1,2] = 0.0

# the faces (triangles)
faces = B['triangles']
# Create meshes
bottom_msh = mesh.Mesh(numpy.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        bottom_msh.vectors[i][j] = vrtb[f[j],:]
        
top_msh = mesh.Mesh(numpy.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        top_msh.vectors[i][j] = vrtt[f[j],:]
# Write meshes to files
bottom_msh.save('bottom_mesh.stl')
top_msh.save('top_mesh.stl')
