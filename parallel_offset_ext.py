#!/usr/bin/python 

import shapefile, sys, fiona, csv
from shapely.geometry import Point, LineString
import matplotlib.pyplot as plt

# This script generates a set of polylines  parallel to all input
# polylines found in the input ESRI shapefile. The lines are generated 
# at the distance offDist from the input line, inside or outside of the 
# input polygons. The output is in  ESRI polyline Z format where Z is 
# (constant) depth. In addition depth is also saved in the record #1 "Depth".

# Required input:
#--------------------------------------------------------------
offDist = 15.0
#----------------
Depth = 0.2
#----------------
WorkDir = "/home/svassili/OpiniconShp/OpiniconModelBuild/"
#InputFile = "Initial_Data/opinicon_perimeter_ed.shp"
InputFile = "TestData/opinicon_perimeter_ed_simpl_0.4-2.shp"
#----------------
# Offset the following shapes to the right, default is offset left.
ShapesRight = [0,12,13,25,31,32,33,36,38,39,40,41,42,45,46,47,50,51,52,53,54]
#----------------
OutFile = "TestData/Perim_ed_simpl_0.4__10m_ext_offset."
#----------------
# If parallel offset algorithm crashes try to adjust coordinate shifts
dX = -200
dY = -200
#---------------------------------------------------------------
InputFile = WorkDir + InputFile
OutFile = WorkDir + OutFile

print "Reading shapefile"
# Read shoreline points
sh = shapefile.Reader(InputFile)
# Shape records
sh_records = sh.shapeRecords()
# Shapes:
sh_shapes = sh.shapes()

# Make temporary coordinates (perimeter + islands) and
# an array of pointers [ind] to starting points of the records
xt=[];yt=[];zt=[];
nShapes=range(0,len(sh_shapes))

j=0; ii=0; ind=[0]
shiftX=sh_shapes[0].points[0][0]
shiftY=sh_shapes[0].points[0][1]

for i in nShapes:
   n=len(sh_shapes[j].points)
   ii=ii+n
   ind.append(ii)
   j+=1
   for k in range(0,n):
       xt.append(sh_shapes[i].points[k][0] - shiftX + dX)
       yt.append(sh_shapes[i].points[k][1] - shiftY + dY) 
       zt.append(0.0)

# Perimeter parallel offset line
#----------------------------------------------------------------------------
# all shapes, left offset is 1, right is 0
print "Generating parallel offset lines"
Shapes=[]
for i in range(0,len(sh_records)):
  Shapes.append(1)
for i in ShapesRight:
  Shapes[i]=0

# Initialize output shapefile
ShapeType=shapefile.POLYLINEZ
w=shapefile.Writer(ShapeType)
w.autobalance=1
w.field("Depth", "F",10,5)

# Compute offset for the following records
#-----------------------------------------
for j in range(1):
#-----------------------------------------
  xf=[]; yf=[]; p=[]
# prepare shapely polyline from the list of coordinates
  for i in range(ind[j],ind[j+1]):
       p.append(Point(xt[i],yt[i]))
       xf.append(xt[i]); yf.append(yt[i])
  line=LineString(p)

# generate parallel offset. 1 - interior offset; 0 - exterior offset
  if Shapes[j] == 0:  
    try:
       offsetLine = line.parallel_offset(distance=offDist, side='left', join_style=2, mitre_limit=5.0)
    except ValueError: 
       sys.exit("Bug in parallel offset, try different dX/dY") 
  else:
    try:  
       offsetLine = line.parallel_offset(distance=offDist, side='right', join_style=2, mitre_limit=5.0)
    except ValueError: 
       sys.exit("Bug in parallel offset, try different dX/dY ")
  ml=0;mli=0    
  try:
  # try to process offset line as a list of polylines
    print  "S" + str(j) + ":",  
    for i in range(len(list(offsetLine))):
       if offsetLine[i].length > ml:
          ml = offsetLine[i].length
          mli = i
       vertex=[]
       x1off,y1off=offsetLine[i].xy
       for k in range(0,len(x1off)):
          vertex.append([x1off[k] + shiftX - dX, y1off[k] + shiftY - dY, Depth])   
       print len(x1off),
       w.record(Depth)
       w.poly([vertex])  
       plt.plot(x1off,y1off,c='r')
       vertex=[]
       x1off,y1off=offsetLine[mli].xy
       for k in range(len(x1off)):
          vertex.append([x1off[k] + shiftX - dX, y1off[k] + shiftY - dY, Depth])   

  except TypeError:
    vertex=[] 
  # if multiline reader fails process offsetLine as a single line 
    x1off,y1off=offsetLine.xy
    print len(x1off),
    for k in range(0, len(x1off)):
       vertex.append([x1off[k] + shiftX - dX, y1off[k] + shiftY - dY, Depth])
    w.record(Depth)
    w.poly([vertex])  
    plt.plot(x1off,y1off,c='r')
   
  plt.plot(xf,yf,c='lightgrey')
  print ""  
  sys.stdout.flush()
   
# <<<<<<< Write out vertices  >>>>>>>
with open('vertices_ext.csv', 'wb') as f:
    writer = csv.writer(f)
    writer.writerows(vertex)
  
for s in w.shapes():
  s.shapeType = ShapeType
w.save(OutFile)

# write projection
with fiona.open(InputFile) as fp:
 prj=open(OutFile+"prj","w")
 prj.write(fp.crs_wkt)
 prj.close()
  
# plot results
plt.show()


