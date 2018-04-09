def align_to_principal(coord,d):
  # xyz - coordinates;
  # d - the number of dimensions for computation of transformation matrix, can be 2 or 3  
  # return numpy array dim[n,3]
  import numpy as np
  print "<<< Computing standard orientation >>>"
  N = coord.shape[0];
  z=[]
    
  # ----- <<< prepare 2D coordinate array >>> ----- #
  if d == 2:
    for i in range(N):
      z.append(coord[i,2])  
      coord[i,2] = 0.
  # ----- <<< compute geometric center >>> ----- #
  center = np.mean(coord, 0)
  print "Geometric center:\n   ", center
  # ----- <<< center coordinates and compute principal axes >>> ----- #
  coord = coord - center
  inertia = np.dot(coord.transpose(), coord)
  e_values, e_vectors = np.linalg.eig(inertia)
  # ----- <<< order eigenvalues (and eigenvectors) >>> ----- #
  order = np.argsort(e_values) 
  eval3, eval2, eval1 = e_values[order]
  axis3, axis2, axis1 = e_vectors[:, order].transpose()
  print "Inertia axes:\n   ", axis1, "\n   ", axis2, "\n   ", axis3
  # ----- <<< orient coordinates >>> ----- #
  Rot = np.array([axis1,axis2,axis3]) # rotation matrix
  for i in range(N):
     if d == 2: # restore z 
          coord[i,2] = z[i] 
     coord[i] = np.dot(Rot,coord[i]) # rotate
  # ----- <<< print stats >>> ----- # 
  print "XYZ axis limits:\n   " , coord[:,0].min(), "to", coord[:,0].max() 
  print "   " , coord[:,1].min(), "to", coord[:,1].max() 
  print "   " , coord[:,2].min(), "to", coord[:,2].max()
  return(coord)

def principal_axes(coord,d):
  # coord - coordinates [nx3];
  # d - dimension for computation of transformation matrix, 2 or 3  
  # returns center [1,3] and rotation matrix [3,3]
  import numpy as np
  print "<<< Computing center and rotation matrix >>>"
  N = coord.shape[0]; A = np.copy(coord)
  # ----- <<< make 3D coordinate array >>> ----- #
  if d == 2:
    for i in range(N): 
      A[i,2] = 0.
  # ----- <<< compute geometric center >>> ----- #
  center = np.mean(A, 0)
  print "... Geometric center:\n   ", center
  # ------- <<< compute principal axes >>> ------ #
  A -= center
  inertia = np.dot(A.transpose(), A)
  e_values, e_vectors = np.linalg.eig(inertia)
  # ----- <<< order eigenvalues (and eigenvectors) >>> ----- #
  order = np.argsort(e_values) 
  eval3, eval2, eval1 = e_values[order]
  axis3, axis2, axis1 = e_vectors[:, order].transpose()
  print "... Inertia axes:\n   ", axis1, "\n   ", axis2, "\n   ", axis3
  # ----- <<< rotation matrix >>> ----- #
  Rot = np.array([axis1,axis2,axis3]) 
  return(center,Rot)

def rotate_coord(coord, center, Rot):
  # coord - coordinates [n,2] or [n,3]
  # center - center of rotation [1,3] 
  # Rot - rotation matrix [3,3]    
  # returns array same dimension as input
  import numpy as np
  print "<<< Transforming coordinates >>>"
  N = coord.shape[0]; n = coord.shape[1]
  A = np.copy(coord)
  if n == 2:
      Rot = Rot[0:2,0:2]
      center = center[0:2]
  # ----- <<< orient coordinates >>> ----- #
  for i in range(N):
       A[i] -= center
       A[i] = np.dot(Rot,A[i]) # rotate
  # ----- <<< print stats >>> ----- # 
  print "... XYZ axis limits:\n   " , A[:,0].min(), "to", A[:,0].max() 
  print "   " , A[:,1].min(), "to", A[:,1].max()
  if n == 3:
      print "   " , A[:,2].min(), "to", A[:,2].max()
  return(A)

