#!/usr/bin/python	

from mergepoints import merge_points
import shapefile, fiona, csv, math, numpy

# This script reads depth points from .csv and shape files, and 
# combines them in a single shapefile. All files should be in
# the same projection. 
#-------------------------------------------------------------------
# Required input:
WorkDir="/home/svassili/OpiniconShp/OpiniconModelBuild/"
PerimeterFile = "Refined_Data/opinicon_perim_&_offset.dbf"
SounderFile1 = "Initial_Data/opinicon_raw_gps_no_duplicates.shp"
csvFile1 = "Initial_Data/20171021_location.csv"
OutFile = "opinicon_combined_bathymetry."
# Fill gaps between points separated by distance bigger than: 
space = 5.0;
# Multiply z b:
zmult = 20.0
# Merge points separated by less than: 
r = 1.0
# Maximum allowed spike 
maxh = 4.0
#------------------------------------------------------------
SounderFile1 = WorkDir + SounderFile1
csvFile1 = WorkDir + csvFile1
PerimeterFile = WorkDir + PerimeterFile
x=[]; y=[]; z=[]; depth=[];


# Read perimeter file
#------------------------------------------------------------
file=PerimeterFile
print "Input:", file
sh = shapefile.Reader(file)
sh_records = sh.shapeRecords()
sh_shapes = sh.shapes()

nShapes=range(len(sh_shapes))
ii=0; j=0; ind=[]; ind.append(0)
segments=[]
holes=[]
count=0
for i in nShapes:
   n=len(sh_shapes[i].points)
   ii=ii+n; hx=0; hy=0
   ind.append(ii)
   for j in range(n):
       tx=sh_shapes[i].points[j][0]
       ty=sh_shapes[i].points[j][1] 
       x.append(tx); y.append(ty)
       hx += tx; hy += ty
   hx /= n; hy /= n
   holes.append([hx,hy])
   z += sh_shapes[i].z
print "Perimeter:", len(x), "depth points\n"

# Build array of segments for constrained Delaunay
for i in range(57):
   for j in range(ind[i],ind[i+1]-1):
      segments.append([j,j+1]) 
   segments.append([j+1,ind[i]])   
   
# Read raw sounder data: 
#------------------------------------------------------------
file=SounderFile1
print "Input:", file
bt = shapefile.Reader(file)
bt_records = bt.shapeRecords()
bt_shapes = bt.shapes()
# depth is record #3 in the "opinicon_raw_gps"
xt = []; yt = []; zt = []
x2 = []; y2 = []; z2 = []
for i in range(len(bt_shapes)):
   for k in range(len(bt_shapes[i].points)):       
       xt.append(bt_shapes[i].points[k][0])
       yt.append(bt_shapes[i].points[k][1])
       zt.append(bt_records[i].record[3])
print "Sounder data:",len(xt),"depth points\n"

# Find and delete outliers
print "... Deleting spikes ..."
print '{0:6} {1:8} {2:6}'.format("     #","     ID","Height")
c=0
for i in range(1,len(xt)-1):
    h = min(zt[i]-zt[i-1], zt[i]-zt[i+1])
    if h > maxh:
        c+=1
        print '{0:6} {1:8} {2:6}'.format(c,i,h)
    else:
        x2.append(xt[i]); y2.append(yt[i]); z2.append(zt[i])

# Merge closely spaced data points
print "... Merging closely spaced data points ..."
while True:
  try:  
    x2,y2,z2 = merge_points(x2,y2,z2,r)
  except TypeError:
    break
x += x2[:]; y += y2[:]; z += z2[:]
print "\n"

# Read depth measurements from csv file
#--------------------------------------------------------
file = csvFile1
print "Input:", file
try:
  with open(file) as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    next(csvReader, None)
    for row in csvReader:
      x.append(float(row[3]))
      y.append(float(row[4]))
      z.append(-float(row[8]))
  print "CSV Table,", int(csvReader.line_num)-1,"depth points\n"
except IOError:
  pass

# Triangulate Planar Straight Line Graph
#----------------------------------------------------------
import triangle
import triangle.plot

ns=len(x)       
XY_S = numpy.ndarray(shape = (ns,2), dtype = float)
for i in range(ns):
   XY_S[i,0] = x[i]; XY_S[i,1] = y[i]; 

ns=len(segments)
SEGM = numpy.ndarray(shape = (ns,2), dtype = int)
for i in range(ns):
   SEGM[i,0] = segments[i][0]; SEGM[i,1] = segments[i][1]; 

ns=56
holes.pop(0)
HOLES = numpy.ndarray(shape = (ns,2), dtype = float)
for i in range(ns):
   HOLES[i,0] = holes[i][0]; HOLES[i,1] = holes[i][1]; 

print "Constrained conforming Delaunay triangulation of the PSLG"   
print "... this may take a while ...\n"
A = dict(vertices=XY_S, segments=SEGM, holes=HOLES)
B = triangle.triangulate(A,'pq7')

# ToDo: Save segments and holes

# Save STL file
V=B['vertices']
T=B['triangles']

import matplotlib.pyplot as plt
triangle.plot.compare(plt, A, B)
plt.show()

#-------------------------------------------------------
# write output shapefile and projection
#--------------------------------------------------------
# Initialize output shapefile
ShapeType=shapefile.POINTZ
w=shapefile.Writer(ShapeType)
w.autobalance=1
w.field("ID", "F",10,5)

print "Output:"
print OutFile+"shp,", len(x), "points"

# write as points
pcount=0
for i in range(len(x)):
  w.point(x[i],y[i],z[i]*zmult)
  w.record(z[i]*zmult)
  pcount += 1
print pcount, "points written"
       
for s in w.shapes():
  s.shapeType = ShapeType
w.save(OutFile)

# Write projection
with fiona.open(PerimeterFile) as fp:
  prj=open(OutFile+"prj","w")
  prj.write(fp.crs_wkt)
  prj.close()




