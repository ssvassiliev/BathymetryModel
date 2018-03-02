#!/usr/bin/python	

from mergepoints import merge_points
import shapefile, fiona, csv, math, numpy

# This script reads depth points from .csv and shape files, and 
# combines them in a single shapefile. All files should be in
# the same projection. 
#-------------------------------------------------------------------
# Required input:
WorkDir="/home/svassili/OpiniconShp/OpiniconModelBuild/"
InputFile1 = "Initial_Data/opinicon_raw_gps_no_duplicates.shp"
InputFile2 = "Initial_Data/20171021_location.csv"
InputFile3 = "Refined_Data/opinicon_perimeter_ed_simpl_0.4-2.shp"
InputFile4 = "Refined_Data/Perim_ed_simpl_0.4__1m_off_0.5m_depth.shp"
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
  for i in range(0,len(bt_shapes)):
    xt = [];yt = []
    for j in range(0,len(bt_shapes[i].points)):
      tx = bt_shapes[i].points[j][0]
      ty = bt_shapes[i].points[j][1]
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
          nn += 1
  print nn, "points added" 
  return(x,y,z) 

# Read raw sounder data: File 1
#-------------------------------------------------------------------
print "Input:"
bt = shapefile.Reader(InputFile1)
bt_records = bt.shapeRecords()
bt_shapes = bt.shapes()
# depth is record #3 in the "opinicon_raw_gps"
for i in range(0,len(bt_shapes)):
   for k in range(0,len(bt_shapes[i].points)):       
       x.append(bt_shapes[i].points[k][0])
       y.append(bt_shapes[i].points[k][1])
       z.append(bt_records[i].record[3])
print "1: Sounder data:",len(x),"depth points"

# Find and delete outliers
xt = []; yt = []; zt = []
print "... deleting spikes ..."
print '{0:6} {1:8} {2:6}'.format("     #","     ID","Height")
c=0
for i in range(1,len(x)-1):
    h = min(z[i]-z[i-1], z[i]-z[i+1])
    if h > maxh:
        c+=1
        print '{0:6} {1:8} {2:6}'.format(c,i,h)
    else:
        xt.append(x[i]); yt.append(y[i]); zt.append(z[i])
x = xt[:]; y = yt[:]; z = zt[:]        

# Merge closely spaced data points
print "... Merging closely spaced data points ..."
while True:
  try:  
    x,y,z = merge_points(x,y,z,r)
  except TypeError:
    break

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

# Read perimeter shape file: File3
#-----------------------------------------------------------
# This file does not have both z coordinate and Depth field, 
# so we add them. We also add points on a straight line where 
# distance between points is bigger than a predefined value
#-----------------------------------------------------------
try:
  np=len(x)  
  bt = shapefile.Reader(InputFile3)
  bt_shapes = bt.shapes()
  for i in range(0,len(bt_shapes)):
    for j in range(0,len(bt_shapes[i].points)):
      tx=bt_shapes[i].points[j][0]
      ty=bt_shapes[i].points[j][1]
      x.append(tx);y.append(ty); z.append(float(0))
  print "3: Perimeter,", len(x) - np, "points",    
  add_points(bt,x,y,z)
except:
  pass

# Read perimeter parallel offset:  File 4 
#-------------------------------------------------------
try:
  np=len(x)  
  bt = shapefile.Reader(InputFile4)
  bt_shapes = bt.shapes()
  bt_records = bt.shapeRecords()
  for i in range(0,len(bt_shapes)):
    for j in range(0,len(bt_shapes[i].points)):
      tx=bt_shapes[i].points[j][0]
      ty=bt_shapes[i].points[j][1]
      tz=float(bt_records[i].record[0])
      x.append(tx);y.append(ty); z.append(tz)
  print "4: P. offset,", len(x) - np, "points",  
  add_points(bt,x,y,z)
except:
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


