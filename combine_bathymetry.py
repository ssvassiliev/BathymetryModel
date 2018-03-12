#!/usr/bin/python	
from mergepoints import merge_points, sqdistance
from stl import mesh
import shapefile, fiona, csv, math, numpy, os, sys

# This script reads depth points from csv and shape files, and 
# combines them in a single shapefile. All files should be in
# the same projection. 
#-------------------------------------------------------------------
# Required input:
WorkDir="/home/svassili/OpiniconShp/OpiniconModelBuild/"
PerimeterFile = "Refined_Data/opinicon_perim_&_offset.dbf"
SounderFile1 = "Initial_Data/opinicon_raw_gps_no_duplicates.shp"
csvFile1 = "Initial_Data/20171021_location.csv"
OutFile = "opinicon_combined_bathymetry."
# z scaling factor.
# Output z will be multiplied by zmult.
zmult = 40.0
# Minimal allowed points separation.
# Points separated by less than r will be merged. 
r = 1.0
# Maximum iterrations
ntries = 20
# Moving average distance 
n_av = 20.0
# Maximum allowed depth peak height.
# Positive and negative peaks with h < maxh and width = 1
# will be replaced by average of z[i-1] and z[i+1]
maxh = 1.0
print "\nMinimal allowed points separation =", r
print "Moving average distance =", n_av
print "Maximum allowed peak height =",maxh
print "Depth scaling factor =",zmult
#------------------------------------------------------------
SounderFile1 = WorkDir + SounderFile1
csvFile1 = WorkDir + csvFile1
PerimeterFile = WorkDir + PerimeterFile
x=[]; y=[]; z=[]; depth=[];
#------------------------------------------------------------
# <<<<<<<<<<<<<<<<< Read perimeter file >>.>>>>>>>>>>>>>>>>>>
#------------------------------------------------------------
file=PerimeterFile
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
print "\n<<< Reading perimeter >>>\n", os.path.basename(file)+",", len(x), "points\n"  
#------------------------------------------------------------
# <<<<<<<<<<<<<<<< Read raw sounder data: >>>>>>>>>>>>>>>>>> 
#------------------------------------------------------------
file=SounderFile1
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
print "<<< Reading sounder data >>>\n", os.path.basename(file)+",", len(xt),"points\n" 
# Find and delete outliers
print "<<< Deleting spikes >>>"
print '{0:6} {1:8} {2:6}'.format("     #","     ID","Height+ " "Height-")
c=0
for i in range(1,len(xt)-1):
   # spikes
    h = min(zt[i]-zt[i-1], zt[i]-zt[i+1])
   # dips
    g = min(-zt[i]+zt[i-1], -zt[i]+zt[i+1])   
    if h > maxh or g > maxh:
        c+=1
        print '{0:6} {1:8} {2:7} {3:7}'.format(c,i,h,g)
        x2.append(xt[i]); y2.append(yt[i]);
        z2.append((zt[i+1]+zt[i-1])*0.5)    
    else:
        x2.append(xt[i]); y2.append(yt[i]); z2.append(zt[i])

# Moving average over a distance n_av
print "\n<<< Computing moving average >>>" 
n_av *= n_av
for i in range(len(x2)):
   j=1;k=0
   tm=0.0
   while i+j < len(x2): 
      if sqdistance(x2,y2,i,i+j) < n_av:
         tm += z2[i+j]
         j += 1
      else:
         if k == i+j:
           break     
         z2[i] += tm
         z2[i] /= j
         k=i+j
         break
   if i + j == len(x):
      break
      
# Merge closely spaced data points
print "\n<<< Merging closely spaced data points >>>"
count=0
while True:
  try:
    if count > ntries:
       print "Error: maximum iterrations reached"
       print "Try to decrease minimal allowed points separation"
       sys.exit()
    x2,y2,z2 = merge_points(x2,y2,z2,r)
    count += 1
  except TypeError:
    break
x += x2[:]; y += y2[:]; z += z2[:]
#--------------------------------------------------------
# <<<<<<< Read depth measurements from csv file >>>>>>>>>
#--------------------------------------------------------
file = csvFile1
try:
  with open(file) as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    next(csvReader, None)
    for row in csvReader:
      x.append(float(row[3]))
      y.append(float(row[4]))
      z.append(-float(row[8]))
  print "\n<<< Reading CSV table>>>\n", os.path.basename(file)+",",\
   int(csvReader.line_num)-1,"depth points\n"
except IOError:
  pass
#---------------------------------------------------------
# <<<<<< Triangulate Planar Straight Line Graph >>>>>>>>>
#---------------------------------------------------------
import triangle
import triangle.plot
print "<<< Constrained conforming Delaunay triangulation >>>"   
print "... this may take a while ...\n"
# Build array of segments
nPoly=len(nShapes)/2
for i in range(nPoly):
   for j in range(ind[i],ind[i+1]-1):
      segments.append([j,j+1]) 
   segments.append([j+1,ind[i]])   
# vertices
ns=len(x)       
XY_S = numpy.ndarray(shape = (ns,2), dtype = float)
for i in range(ns):
   XY_S[i,0] = x[i]; XY_S[i,1] = y[i]; 
#segments
ns=len(segments)
SEGM = numpy.ndarray(shape = (ns,2), dtype = int)
for i in range(ns):
   SEGM[i,0] = segments[i][0]; SEGM[i,1] = segments[i][1]; 
# number of holes (islands)
ns=nPoly-1
# perimeter comes first, so we skip it
holes.pop(0)
HOLES = numpy.ndarray(shape = (ns,2), dtype = float)
for i in range(ns):
   HOLES[i,0] = holes[i][0]; HOLES[i,1] = holes[i][1]; 

A = dict(vertices=XY_S, segments=SEGM, holes=HOLES)
B = triangle.triangulate(A,'pq0')
#----------------------------------------------------
# <<<<<<<<<<<<<<<<< Write stl files >>>>>>>>>>>>>>>>>
#----------------------------------------------------
ns=len(x)
na = len(B['vertices'])
# the existing vertices
vrtb = numpy.ndarray(shape = (na,3), dtype = float)
vrtt = numpy.ndarray(shape = (na,3), dtype = float)
for i in range(ns):
   vrtb[i,0] = x[i]; vrtb[i,1] = y[i]; vrtb[i,2] = z[i]*zmult;
   vrtt[i,0] = x[i]; vrtt[i,1] = y[i]; vrtt[i,2] = 0.0;
# the new vertices
for i in range(ns,na):
   vrtb[i,0] = vrtt[i,0] = B['vertices'][i][0];   
   vrtb[i,1] = vrtt[i,1] = B['vertices'][i][1];
   vrtb[i,2] = vrtt[i,2] = 0.0;
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
#----------------------------------------------------------
# <<<<<<<< Write output shapefile and projection >>>>>>>>>>
#----------------------------------------------------------
# Initialize output shapefile
ShapeType=shapefile.POINTZ
w=shapefile.Writer(ShapeType)
w.autobalance=1
w.field("ID", "F",10,5)
print "<<< Output >>>"
print OutFile+"shp,", len(x), "points"
# write as points
pcount=0
for i in range(len(x)):
  w.point(x[i],y[i],z[i]*zmult)
  w.record(z[i]*zmult)
  pcount += 1
for s in w.shapes():
  s.shapeType = ShapeType
w.save(OutFile)
# Write projection
with fiona.open(PerimeterFile) as fp:
  prj=open(OutFile+"prj","w")
  prj.write(fp.crs_wkt)
  prj.close()

#import matplotlib.pyplot as plt
#triangle.plot.compare(plt, A, B)
#plt.show()


