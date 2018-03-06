#!/usr/bin/python	

from mergepoints import merge_points
import shapefile, fiona, csv, math, numpy

# This script reads depth points from .csv and shape files, and 
# combines them in a single shapefile. All files should be in
# the same projection. 
#-------------------------------------------------------------------
# Required input:
WorkDir="/home/svassili/OpiniconShp/OpiniconModelBuild/"
InputFile1 = "Refined_Data/opinicon_perimeter_ed_simpl_0.4-2.shp"
InputFile2 = "Refined_Data/Perim_ed_simpl_0.4__1m_off_0.5m_depth.shp"
InputFile3 = "Initial_Data/opinicon_raw_gps_no_duplicates.shp"
InputFile4 = "Initial_Data/20171021_location.csv"
OutFile = "opinicon_combined_bathymetry."
# Fill gaps between points separated by distance bigger than: 
space = 5.0;
# Multiply z by:
zmult = 20.0
# Merge points separated by less than: 
r = 1.0
# Maximum allowed spike 
maxh = 4.0
#-------------------------------------------------------------------
InputFile1 = WorkDir + InputFile1
InputFile2 = WorkDir + InputFile2
InputFile3 = WorkDir + InputFile3
InputFile4 = WorkDir + InputFile4
x=[]; y=[]; z=[]; depth=[];

def add_points(bt,x,y,z):
  np=len(x); nn=0  
  bt_shapes = bt.shapes()
  for i in range(len(bt_shapes)):
    xt = []; yt = []; zt=[]
    for j in range(len(bt_shapes[i].points)):
      tx = bt_shapes[i].points[j][0]
      ty = bt_shapes[i].points[j][1]
      tz=float(bt_records[i].record[0])
      xt.append(tx); yt.append(ty); zt.append(tz);
    for k in range(len(xt)-1):
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
        for l in range(n):
          x.append(xt[k]+dx*l*sp)
          y.append(yt[k]+dy*l*sp)
          z.append(zt[k])
          nn += 1
  print nn, "points added" 
  return(x,y,z) 

def add_points_2(bt,x,y,z,ii,ind):
  np=len(x); nn=0  
  bt_shapes = bt.shapes()
  for i in range(len(bt_shapes)):
    xt = []; yt = []; zt=[]
    ns=len(bt_shapes[i].points)
    for j in range(ns-1):
      tx1 = bt_shapes[i].points[j][0]
      ty1 = bt_shapes[i].points[j][1]
      tz1 = float(bt_records[i].record[0])
      tx2 = bt_shapes[i].points[j+1][0]
      ty2 = bt_shapes[i].points[j+1][1]
      x.append(tx1); y.append(ty1); z.append(tz1);ii+=1;
      if j == ns-1:
        dx=x[0]-x[j]; dy=y[0]-y[j]
      else:
        dx=x[j+1]-x[j]; dy=y[j+1]-y[j]
      dist=math.sqrt(dx*dx+dy*dy)
      idist=1/dist
      dx=dx*idist; dy=dy*idist
      if dist > space:
        n=int(math.ceil(dist/space))
        sp=dist/n
        for l in range(n):
          x.append(xt[j]+dx*l*sp)
          y.append(yt[j]+dy*l*sp)
          z.append(zt[j])
          ii += 1
          nn += 1
  ind.append(ii)
  print nn, "points added" 
  return(x,y,z) 




# Read perimeter shape file: File1
#-----------------------------------------------------------
# This file does not have both z coordinate and Depth field, 
# so we add them. We also add points on a straight line where 
# distance between points is bigger than a predefined value
#-----------------------------------------------------------


ii=0; ind=[]
try: 
  bt = shapefile.Reader(InputFile1)
  bt_shapes = bt.shapes()
  bt_records = bt.shapeRecords()
  for i in range(len(bt_shapes)):
    n=len(bt_shapes[i].points)
    ii += n
    ind.append(ii)
    for j in range(n):
      tx=bt_shapes[i].points[j][0]
      ty=bt_shapes[i].points[j][1]
      tz=float(bt_records[i].record[0])
      x.append(tx);y.append(ty); z.append(tz)
  print "1: Perimeter,", len(x), "points",    
  add_points(bt,x,y,z)
except:
  pass

# Read perimeter parallel offset:  File 2 
#-------------------------------------------------------
try:
  np=len(x)  
  bt = shapefile.Reader(InputFile2)
  bt_shapes = bt.shapes()
  bt_records = bt.shapeRecords()
  for i in range(len(bt_shapes)):
    n=len(bt_shapes[i].points)
    ii += n
    ind.append(ii)
    for j in range(len(bt_shapes[i].points)):
      tx=bt_shapes[i].points[j][0]
      ty=bt_shapes[i].points[j][1]
      tz=float(bt_records[i].record[0])
      x.append(tx);y.append(ty); z.append(tz)
  print "2: P. offset,", len(x) - np, "points",  
  add_points(bt,x,y,z)
except:
  pass

# Read raw sounder data: File 3
#----------------------------------------------------------
print "Input:"
bt = shapefile.Reader(InputFile3)
bt_records = bt.shapeRecords()
bt_shapes = bt.shapes()
# depth is record #3 in the "opinicon_raw_gps"
xt = []; yt = []; zt = []
x2 = []; y2 = []; z2 = []
for i in range(len(bt_shapes)):
   for k in range(0,len(bt_shapes[i].points)):       
       xt.append(bt_shapes[i].points[k][0])
       yt.append(bt_shapes[i].points[k][1])
       zt.append(bt_records[i].record[3])
print "3: Sounder data:",len(xt),"depth points"

# Find and delete outliers
print "... deleting spikes ..."
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
ii += len(x2)
ind.append(ii)
  
# Read depth measurements from csv file: File 4
#--------------------------------------------------------
try:
  with open(InputFile4) as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    next(csvReader, None)
    for row in csvReader:
      ii += 1
      x.append(float(row[3]))
      y.append(float(row[4]))
      z.append(-float(row[8]))
  print "4: CSV Table,", int(csvReader.line_num)-1,"depth points"
  ind.append(ii)
except IOError:
  pass

# write output shapefile and projection
#--------------------------------------------------------
# Initialize output shapefile
ShapeType=shapefile.POINTZ
w=shapefile.Writer(ShapeType)
w.autobalance=1
w.field("Depth", "F",10,5)

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


