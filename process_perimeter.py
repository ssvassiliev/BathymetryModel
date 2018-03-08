#!/usr/bin/python	

import shapefile, fiona, csv, math, numpy

# This script reads depth points from perimeter and perimeter
# parallel offset files, and combines them in a single shapefile
# filling gaps in straight line segments. Both files should be in
# the same projection. 
#-------------------------------------------------------------------
# Required input:
WorkDir="/home/svassili/OpiniconShp/OpiniconModelBuild/"
InputFile1 = "Refined_Data/opinicon_perimeter_ed_simpl_0.4-2.shp"
InputFile2 = "Refined_Data/Perim_ed_simpl_0.4__1m_off_0.5m_depth.shp"
OutFile = "opinicon_perim_&_offset."
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

# add_points function takes data points from shapefile reader object
# bt and appends them to arrays (x,y,z). It adds points on a straight
# line where distance between dat points is bigger than a predefined
# (spacing) parameter. Returns updated number of data points (ii) and
# updated array of indexes (ind) pointing to the last element of each
# read shape. This code works only for shape records in which the last
# point is the same as the first one, and depth is stored in the first
# record.  
#-----------------------------------------------------------------
def add_points(bt,x,y,z,ii,ind,spacing):
  import shapefile, math
  n_read=0; n_added=0; n_prev=ii
  bt_shapes = bt.shapes()
  bt_records = bt.shapeRecords()
  for i in range(len(bt_shapes)):
    ns=len(bt_shapes[i].points)
    tx0=bt_shapes[i].points[0][0]
    ty0=bt_shapes[i].points[0][1]
    for j in range(0,ns-1):     
      tx1 = bt_shapes[i].points[j][0]
      ty1 = bt_shapes[i].points[j][1]
      tz1 = float(bt_records[i].record[0])
      x.append(tx1); y.append(ty1); z.append(tz1); n_read += 1  
      tx2 = bt_shapes[i].points[j+1][0]
      ty2 = bt_shapes[i].points[j+1][1]    
      dx = tx2 - tx1; dy = ty2 - ty1     
      dist=math.sqrt(dx*dx+dy*dy)
      idist=1/dist; dx=dx*idist; dy=dy*idist
      if dist > spacing:
        n=int(math.ceil(dist/spacing))
        sp=dist/n
        for l in range(1,n):
          x.append(tx1+dx*l*sp)
          y.append(ty1+dy*l*sp)
          z.append(z[j])
          n_added += 1
    ind.append(n_prev+n_read+n_added)
  print n_read, "points",  
  print n_added, "points added"
  ii = n_prev + n_read + n_added
  return(x,y,z,ii,ind)

#-------------------------------------------------------
#              Read perimeter
#-------------------------------------------------------
ii=0; ind=[]
ind.append(0)
try:
  print "1: Perimeter,",  
  bt = shapefile.Reader(InputFile1)    
  x,y,z,ii,ind = add_points(bt,x,y,z,ii,ind,space)
except:
  pass
#------------------------------------------------------
#       Read perimeter parallel offset
#------------------------------------------------------
try:
  print "2: Perimeter,",  
  bt = shapefile.Reader(InputFile2)    
  x,y,z,ii,ind = add_points(bt,x,y,z,ii,ind,space)
except:
  pass
#-------------------------------------------------------
# write output shapefile and projection
#--------------------------------------------------------
# Initialize output shapefile
ShapeType=shapefile.POLYLINEZ
w=shapefile.Writer(ShapeType)
w.autobalance=1
w.field("ID", "F",10,5)

print "Output:"
print OutFile+"shp,", len(x), "points"

# write as polyline
vcount=0
for i in range(0,len(ind)-1):
  vertex=[] 
  for j in range(ind[i],ind[i+1]):
    vertex.append([x[j], y[j], z[j]*zmult, 0])
    vcount += 1
  w.record(i)
  w.line([vertex])
print vcount, "points written"
       
for s in w.shapes():
  s.shapeType = ShapeType
w.save(OutFile)

# Write projection
with fiona.open(InputFile1) as fp:
  prj=open(OutFile+"prj","w")
  prj.write(fp.crs_wkt)
  prj.close()


