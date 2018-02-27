#!/usr/bin/python
import shapefile, fiona, csv, math, os

# This script adds points on a straight line where 
# distance between points is bigger than a threshold
# Each point will be saved in a separate record

#-------------------------------------------------------------------
# Required input:
WorkDir="/home/svassili/OpiniconShp/OpiniconModelBuild/"
InputFile = "Refined_Data/opinicon_perimeter_ed_simpl_0.4.shp"
OutFile = os.path.splitext(InputFile)[0]+"_extra_points."
# Fill gaps between points separated by distance bigger than 'space' 
space = 5.0; 
#-------------------------------------------------------------------
InputFile = WorkDir + InputFile

# Initialize output shapefile
ShapeType=shapefile.POINTZ
w=shapefile.Writer(ShapeType)
w.autobalance=1
w.field("Depth", "F",10,5)
x=[]; y=[]; z=[]; depth=[];

# Read shape file
#-----------------------------------------------------------
np=len(x); nn=0
try:
    bt = shapefile.Reader(InputFile)
    bt_records = bt.shapeRecords()
    bt_shapes = bt.shapes()
    try:
       bt_records[1].z
       depth="Z"
       print "Depth will be read from Z coordinate"
    except AttributeError:  
       print "Input file does not have z coordinate."
       print "Depth will be read from the first field of records"
       depth="F1"

    for i in range(0,len(bt_shapes)):
      xt=[];yt=[];zt=[]
      for j in range(0,len(bt_shapes[i].points)):
          tx=bt_shapes[i].points[j][0]
          ty=bt_shapes[i].points[j][1]
          if depth == "Z":
              tz=bt_records[i].z
          if depth == "F1":
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

    print "File 4:",len(x)-np,"depth points, added", nn, "points" 
except:
   pass

# write output shapefile and projection
#--------------------------------------------------------
ShapeType=shapefile.POINTZ
w=shapefile.Writer(ShapeType)
w.autobalance=1
w.field("Depth", "F",10,5)

for i in range(0,len(x)):
   w.point(x[i],y[i],z[i])
   w.record(z[i]) 
for s in w.shapes():
  s.shapeType = ShapeType
w.save(OutFile)

# Write projection
with fiona.open(InputFile) as fp:
 prj=open(OutFile+"prj","w")
 prj.write(fp.crs_wkt)
 prj.close()

