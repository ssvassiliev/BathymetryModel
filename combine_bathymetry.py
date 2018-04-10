#!/usr/bin/python

from pyproj import Proj, transform
import scipy, numpy, pandas
from scipy.spatial.distance import cdist
from mergepoints import merge_points, sqdistance
from orient import principal_axes, rotate_coord
from stl import mesh
import shapefile, fiona, csv, math, numpy, scipy, os, sys
print "\n************* Combine bathymetry **************"
# This script reads depth points from csv and shape files, and 
# combines them in a single shapefile. All files should be in
# the same projection. 
#-------------------------------------------------------------------
# Required input:
WorkDir = os.getcwd()+"/"
PerimeterFile = WorkDir+"Opinicon/Output/opinicon_perim_and_offset-3.shp"
SounderFile1 = WorkDir+"Opinicon/Data/opinicon_raw_gps_no_duplicates.shp"
csvFile1 = WorkDir+"Opinicon/Data/20171021_location.csv"
xlFile1= WorkDir+"Opinicon/Data/L.opinicon-04-2018.xlsx"
OutFile = WorkDir+"Opinicon/Output/opinicon_combined_bathymetry."
verticesFile = WorkDir+'Opinicon/Output/vertices.csv'
vertices_extFile = WorkDir+'Opinicon/Output/vertices_ext.csv'
holesFile = WorkDir+'Opinicon/Output/holes.csv'
segmentsFile = WorkDir+'Opinicon/Output/segments.csv'
#---- z scaling factor, output will be multiplied by zmult.
zmult = 40.0
#----- Minimal allowed points separation, points separated by less than r will be merged. 
r = 1.0
#----- Maximum merging iterrations
ntries = 20
#----- Moving average distance 
n_av = 10.0
#----- Maximum allowed peak height.
# Positive and negative peaks with h < maxh and width = 1
# will be replaced by average of z[i-1] and z[i+1]
maxh = 1.0
#----- Triangulation options: http://dzhelil.info/triangle/
tri = 'pq20'
#----- Inverse distance interpolation:
invp = 2  # power
intn = 10 # number of nearest neighbours 
#----- Flip coordinates for 3D printing
flip = True

print "... Minimal allowed points separation =", r
print "... Moving average distance =", n_av
print "... Maximum allowed peak height =",maxh
print "... Depth scaling factor =",zmult
#------------------------------------------------------------
x=[]; y=[]; z=[]; depth=[];
#------------------------------------------------------------
# <<<<<<<<<<<<<<<<< Read perimeter file >>.>>>>>>>>>>>>>>>>>>
#------------------------------------------------------------
sh = shapefile.Reader(PerimeterFile)
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
print "<<< Reading perimeter >>>\n...", os.path.basename(PerimeterFile)+",", len(x), "points"
print "... min =", min(z), "max = ", max(z)
#------------------------------------------------------------
# <<<<<<<<<<<<<<<< Read raw sounder data: >>>>>>>>>>>>>>>>>> 
#------------------------------------------------------------
bt = shapefile.Reader(SounderFile1)
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

#---------------------------------------
#<<<<< remove suspicious points >>>>>> 
i1=15993;ni=33
for i in range(ni):
  del xt[i1]; del yt[i1]; del zt[i1]
#---------------------------------------

print "<<< Reading sounder data >>>\n...", os.path.basename(SounderFile1)+",", len(xt),"points" 
# Find and delete outliers
print "<<< Removing outliers >>>"
# print '{0:6} {1:8} {2:6}'.format("     #","     ID","Height+ " "Height-")
c=0
for i in range(1,len(xt)-1):
   # spikes
    h = min(zt[i] - zt[i-1], zt[i] - zt[i+1])
   # dips
    g = min(-zt[i] + zt[i-1], -zt[i] + zt[i+1])   
    if h > maxh or g > maxh:
        c+=1
#        print '{0:6} {1:8} {2:7} {3:7}'.format(c,i,h,g)
        x2.append(xt[i]); y2.append(yt[i]);
        z2.append((zt[i+1] + zt[i-1])*0.5)    
    else:
        x2.append(xt[i]); y2.append(yt[i]); z2.append(zt[i])
print "... removed",c,"outliers" 
# Moving average over a distance n_av
print "<<< Computing moving average >>>" 
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
print "<<< Merging closely spaced data points >>>"
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
try:
  with open(csvFile1) as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    next(csvReader, None)
    for row in csvReader:
      x.append(float(row[3]))
      y.append(float(row[4]))
      z.append(-float(row[8]))
  print "<<< Reading CSV table>>>\n...", os.path.basename(csvFile1)+",",\
   int(csvReader.line_num)-1,"depth points"
except IOError:
  pass

#--------------------------------------------------------
# <<<<<<< Read depth measurements from excel file >>>>>>>>>
#--------------------------------------------------------
# Projections:
wgs84 = Proj('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')##4326
epsg26918 = Proj(init='epsg:26918')

try:
    count=0
    wb = pandas.read_excel(xlFile1)
    t0=list(wb[' longitude'])
    t1=list(wb[' latitude'])
    t2=list(wb['water.depth'])
    xt,yt,zt=transform(wgs84,epsg26918,t0,t1,t2,radians=False)
    n=len(xt)
    for i in range(n):
          x.append(float(xt[i]))
          y.append(float(yt[i]))
          z.append(-float(zt[i]))
          count+=1
    print "<<< Reading excel table>>>\n...", os.path.basename(xlFile1)+",",\
      count,"depth points"
except IOError:
    pass
 
#------------------------------------------------------------
#<<<<<<<<<<<<<<< Add bounding box >>>>>>>>>>>>>>>>>>>
#marg=100
#x1=min(x);x2=max(x);y1=min(y);y2=max(y)
#x.append(x1-marg);y.append(y1-marg);z.append(0.0)
#x.append(x1-marg);y.append(y2+marg);z.append(0.0)
#x.append(x2+marg);y.append(y2+marg);z.append(0.0)
#x.append(x2+marg);y.append(y1-marg);z.append(0.0)
#------------------------------------------------------------

#---------------------------------------------------------
# <<<<<< Triangulate Planar Straight Line Graph >>>>>>>>>
#---------------------------------------------------------
import triangle
import triangle.plot
print "<<< Constrained conforming Delaunay triangulation >>>"   
# segments from offset lines
nPoly=len(nShapes)/2
# offset line segments: 
for i in range(len(nShapes)/2, len(nShapes)):
   for j in range(ind[i],ind[i+1]-1):
      segments.append([j,j+1]) 
   segments.append([j+1,ind[i]])
   
ns=len(segments)
SEGM = numpy.ndarray(shape = (ns,2), dtype = int)
for i in range(ns):
   SEGM[i,0] = segments[i][0]; SEGM[i,1] = segments[i][1]; 

# vertices
ns=len(x)       
XY_S = numpy.ndarray(shape = (ns,2), dtype = float)
for i in range(ns):
   XY_S[i,0] = x[i]; XY_S[i,1] = y[i];    
# number of holes (islands)
ns=nPoly-1
# perimeter comes first, so we skip it
holes.pop(0)
HOLES = numpy.ndarray(shape = (ns,2), dtype = float)
for i in range(ns):
   HOLES[i,0] = holes[i][0]; HOLES[i,1] = holes[i][1]; 

if flip == True:
   center = numpy.mean(XY_S, 0)
   Rot = numpy.ndarray(shape = (3,3), dtype = float)
   Rot[0]=[-1,0,0]
   Rot[1]=[ 0,1,0]
   Rot[2]=[ 0,0,1]
   XY_S=rotate_coord(XY_S,center,Rot)
   HOLES=rotate_coord(HOLES,center,Rot)

A = dict(vertices=XY_S, segments=SEGM, holes=HOLES)
B = triangle.triangulate(A,tri)

#----------------------------------------------------
# <<<<<<<<<< Prepare to write stl files >>>>>>>>>>>>
#----------------------------------------------------
ns=len(x)
na = len(B['vertices'])

# allocate arrays for vertices
vrtb = numpy.ndarray(shape = (na,3), dtype = float)
vrtt = numpy.ndarray(shape = (na,3), dtype = float)
rpt = numpy.ndarray(shape = (1,2), dtype = float)

# fill arrays with existing data 
for i in range(ns):
      vrtb[i,0] = XY_S[i,0]; vrtb[i,1] = XY_S[i,1];vrtb[i,2] = z[i]*zmult; # bottom surface
      vrtt[i,0] = XY_S[i,0]; vrtt[i,1] = XY_S[i,1]; vrtt[i,2] = 0.0;  # top surface

# Interpolate depth of the new points using inverse distance algorithm       
print "<<< Inverse distance weighting >>>"
print "... interpolating depth of", na-ns, "vertices ..."
print "... this may take a while ..."

for j in range(ns,na):
   # the new vertices
   vrtb[j,0] = vrtt[j,0] = B['vertices'][j][0];   
   vrtb[j,1] = vrtt[j,1] = B['vertices'][j][1];
   vrtt[j,2] = 0.0;  
   rpt[0][0] = vrtb[j,0]; rpt[0][1] = vrtb[j,1];  
   dist = scipy.spatial.distance.cdist(rpt,XY_S,'sqeuclidean')
   ds = numpy.argsort(dist)
   sm = 0; mu = 0
   for k in range(intn):
      if abs(z[ds[0][k]]) > 1e-8 : # exclude zero lines from IDW
        di = math.sqrt(dist[0][ds[0][k]])
        wu = pow(di,-invp)
        sm += wu
        mu += z[ds[0][k]]*wu
   if sm != 0.0:     
        vrtb[j,2] = mu*zmult/sm
   else:
        vrtb[j,2] = 0.0

if flip == True:
   for j in range(na):
      vrtb[j,2] *= -1

# the faces (triangles)
faces = B['triangles']

# align with principal axes
center,Rot = principal_axes(vrtb,2)
vrtb = rotate_coord(vrtb,center,Rot)
vrtt = rotate_coord(vrtt,center,Rot)  
HOLES = rotate_coord(HOLES,center,Rot)

# <<<<<<<<< Create meshes >>>>>>>>>>
#-------------------------------------------------------------
bottom_msh = mesh.Mesh(numpy.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        bottom_msh.vectors[i][j] = vrtb[f[j],:]
        
top_msh = mesh.Mesh(numpy.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        top_msh.vectors[i][j] = vrtt[f[j],:]
# Write meshes to files
bottom_msh.save('bottom_mesh_.stl')
top_msh.save('top_mesh_.stl')

#----------------------------------------------------
# <<<<<<< Write out vertices, segments, holes >>>>>>>
with open(verticesFile, 'wb') as f:
    writer = csv.writer(f)
    writer.writerows(vrtb)
segments=[]
# Build array of segments from zero lines
for i in range(nPoly):
   for j in range(ind[i],ind[i+1]-1):
      segments.append([j,j+1]) 
   segments.append([j+1,ind[i]])
#---------------------------------   
# <<< Append bounding box >>>
#ns=len(x)
#segments.append([ns-1,ns-2])
#segments.append([ns-2,ns-3])
#segments.append([ns-3,ns-4])
#segments.append([ns-1,ns-4])
#---------------------------------
ns=len(segments)
SEGM2 = numpy.ndarray(shape = (ns,2), dtype = int)
for i in range(ns):
   SEGM2[i,0] = segments[i][0]; SEGM2[i,1] = segments[i][1]; 
    
with open(segmentsFile, 'wb') as f:
    writer = csv.writer(f)
    writer.writerows(SEGM2)
with open(holesFile, 'wb') as f:
    writer = csv.writer(f)
    writer.writerows(HOLES)
    
    
#----------------------------------------------------------
# <<<<<<<< Write output shapefile and projection >>>>>>>>>>
#----------------------------------------------------------
# Initialize output shapefile
ShapeType=shapefile.POINTZ
w=shapefile.Writer(ShapeType)
w.autobalance=1
w.field("ID", "F",10,5)
print "<<< Output >>>\n...", os.path.basename(OutFile)+"shp,", len(x), "points"
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


