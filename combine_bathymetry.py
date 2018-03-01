#!/usr/bin/python	

import shapefile, fiona, csv, math, numpy, scipy, timeit
from scipy.spatial.distance import pdist, squareform

# This script reads depth points from .csv and shape files, and 
# combines them in a single shapefile. All files should be in
# the same projection. 
#-------------------------------------------------------------------
# Required input:
WorkDir="/home/svassili/OpiniconShp/OpiniconModelBuild/"
InputFile1 = "Initial_Data/opinicon_raw_gps_no_duplicates.shp"
InputFile2 = "Initial_Data/20171021_location.csv"
InputPerim = "Refined_Data/opinicon_perimeter_ed_simpl_0.4-2.shp"
InputFile4 = "Refined_Data/Perim_ed_simpl_0.4__1m_off_0.5m_depth.shp"
OutFile = "opinicon_combined_bathymetry."
# Fill gaps between points separated by distance bigger than 'space'
space = 5.0;
# Multiply z by
zmult=20.0
#-------------------------------------------------------------------
InputFile1 = WorkDir + InputFile1
InputFile2 = WorkDir + InputFile2
InputPerim = WorkDir + InputPerim
InputFile4 = WorkDir + InputFile4

print "Input:"
# Initialize output shapefile
ShapeType=shapefile.POINTZ
w=shapefile.Writer(ShapeType)
w.autobalance=1
w.field("Depth", "F",10,5)
x=[]; y=[]; z=[]; depth=[];

def sq2cond(i, j, n):
    assert i != j, "no diagonal elements in condensed matrix"
    if i < j:
        i, j = j, i
    return n*j - j*(j+1)/2 + i - 1 - j

def row_col_from_condensed_index(d,i):
    b = 1 - 2*d 
    x = math.floor((-b - math.sqrt(b*b - 8*i))*0.5)
    g = i + x*(b + x + 2)*0.5 + 1 
    return (x,g)  


def vec_row_col(d,i):                                                               
  b = 1 - 2* d
  x = (numpy.floor((-b - numpy.sqrt(b*b - 8*i))*0.5)).astype(int)
  g = (i + x*(b + x + 2)*0.5 + 1).astype(int)
  if i.shape:                                                                     
    return zip(x,g)                                                             
  else:                                                                           
    return (x,g) 

def sqdistance(p1,p2):
  import math
  sqd = (x[p1]-x[p2])*(x[p1]-x[p2]) + (y[p1]-y[p2])*(y[p1]-y[p2]) 
  return sqd



# Read raw sounder data: File 1
#--------------------------------------------------------
bt = shapefile.Reader(InputFile1)
bt_records = bt.shapeRecords()
bt_shapes = bt.shapes()


# depth is record #3 in the "opinicon_raw_gps"
for i in range(0,len(bt_shapes)):
   for k in range(0,len(bt_shapes[i].points)):       
       x.append(bt_shapes[i].points[k][0])
       y.append(bt_shapes[i].points[k][1])
       z.append(bt_records[i].record[3])

print "1: Sounder data:",len(x),"depth points,",

# Compute distance matrix 
n = len(x); r = 2.0; 
XY = numpy.ndarray(shape = (n,2), dtype = float)
for i in range(0,n):
  XY[i,0] = x[i]
  XY[i,1] = y[i]
D = scipy.spatial.distance.pdist(XY, 'sqeuclidean')
ix = numpy.where(D < r*r)
rc = vec_row_col(n,ix[0])
print len(rc), "groups with separation <", r, "m"

# Make a list of point groups 
cl=[];groups=[]
pairs=list(rc)
cl.append(pairs[0])
for i in range(1,len(pairs)):
    if pairs[i][0] == pairs[i-1][0]:
       cl.append(pairs[i])
    else:
       groups.append(cl) 
       cl=[]
       cl.append(pairs[i])

# Compute centroids and mark point for deletion
cent_x=[]; cent_y=[]; cent_z=[]; g_i=[True]*n; i_done=[False]*n;
xx=[]; yy=[]; zz=[]
for i in range(0,len(groups)):
    ii=groups[i][0][0]
    if i_done[ii]:
      continue
    cx = x[ii]; cy = y[ii]; cz = z[ii]
    g_i[ii]=False; i_done[ii]=True
    for j in range(0,len(groups[i])):
      ij=groups[i][j][1]
      if i_done[ij]:
        g_i[ii]=True  
        continue
      cx += x[ij]; cy += y[ij]; cz += z[ij]
      g_i[ij]=False; i_done[ij]=True
    div=(len(groups[i])+1) 
    cent_x.append(cx/div); cent_y.append(cy/div); cent_z.append(cz/div)
# Delete groups and replace them with centroids
for i in range(0,n):
   if g_i[i]:
       xx.append(x[i])
       yy.append(y[i])
       zz.append(z[i])
for i in range(0,len(cent_x)):
    xx.append(cent_x[i])
    yy.append(cent_y[i])
    zz.append(cent_z[i])
x=xx[:]; y=yy[:]; z=zz[:]

# Read depth measurements from csv file: File 2
#--------------------------------------------------------
try:
  with open(InputFile2) as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    next(csvReader, None)
    for row in csvReader:
      x.append(float(row[3]))
      y.append(float(row[4]))
      z.append(-float(row[8]))
  print "2: CSV Table,", int(csvReader.line_num)-1,"depth points"
except IOError:
  pass

# Read perimeter shape file: InputPerim
#-----------------------------------------------------------
# This file does not have both z coordinate and Depth field, 
# so we add them. We also add points on a straight line where 
# distance between points is bigger than a predefined value
#-----------------------------------------------------------

np=len(x); nn=0
try:
  bt = shapefile.Reader(InputPerim)
  bt_records = bt.shapeRecords()
  bt_shapes = bt.shapes()
  for i in range(0,len(bt_shapes)):
    xt=[];yt=[]
    for j in range(0,len(bt_shapes[i].points)):
      tx=bt_shapes[i].points[j][0]
      ty=bt_shapes[i].points[j][1]
      x.append(tx);y.append(ty); z.append(float(0.0))
      xt.append(tx); yt.append(ty);
    for k in range(0,len(xt)-1):
      if k == len(xt)-1:
        dx=xt[0]-xt[k]; dy=yt[0]-yt[k]
      else:
        dx=xt[k+1]-xt[k]; dy=yt[k+1]-yt[k]
      dist=math.sqrt(dx*dx+dy*dy)
      idist=1/dist
      dx=dx*idist; dy=dy*idist
      if dist > space:
        n=int(math.ceil(dist/space))
        sp=dist/n
        for l in range(0,n):
          x.append(xt[k]+dx*l*sp)
          y.append(yt[k]+dy*l*sp)
          z.append(float(0.0))
          nn=nn+1
  print "3: Shoreline,",len(x)-np,"points,", nn, "points added" 
except:
  pass

# Read perimeter parallel offset:  File 4 
#-------------------------------------------------------
np=len(x); nn=0
try:
  bt = shapefile.Reader(InputFile4)
  bt_records = bt.shapeRecords()
  bt_shapes = bt.shapes()
  for i in range(0,len(bt_shapes)):
    xt=[];yt=[];zt=[]
    for j in range(0,len(bt_shapes[i].points)):
      tx=bt_shapes[i].points[j][0]
      ty=bt_shapes[i].points[j][1]
      tz=float(bt_records[i].record[0])
      x.append(tx); y.append(ty); z.append(tz)
      xt.append(tx); yt.append(ty); zt.append(tz)
    for k in range(0,len(xt)-1):
      if k == len(xt)-1:
        dx=xt[0]-xt[k]; dy=yt[0]-yt[k]
      else:
        dx=xt[k+1]-xt[k]; dy=yt[k+1]-yt[k]
      dist=math.sqrt(dx*dx+dy*dy)
      idist=1/dist
      dx=dx*idist; dy=dy*idist
      if dist > space:
        n=int(math.ceil(dist/space))
        sp=dist/n
        for l in range(0,n):
           x.append(xt[k]+dx*l*sp)
           y.append(yt[k]+dy*l*sp)
           z.append(zt[k])
           nn=nn+1
  print "4: Shoreline offset,",len(x)-np,"points,", nn, "points added" 
except:
  pass

# write output shapefile and projection
#--------------------------------------------------------
ShapeType=shapefile.POINTZ
w=shapefile.Writer(ShapeType)
w.autobalance=1
w.field("Depth", "F",10,5)

print "Output:"
print OutFile+"shp,", len(x), "points" 
for i in range(0,len(x)):
  w.point(x[i],y[i],z[i]*zmult)
  w.record(z[i]*zmult)
for s in w.shapes():
  s.shapeType = ShapeType
w.save(OutFile)

# Write projection
with fiona.open(InputFile1) as fp:
  prj=open(OutFile+"prj","w")
  prj.write(fp.crs_wkt)
  prj.close()
