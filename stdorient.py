def align_to_principal(coord,d):
  # xyz - coordinates;  d (2 or 3) - the number of dimensions used
  # to compute transformation matrix, transformation itself is always 3D
  # return numpy array dim[n,3]
  import numpy as np
  print "<<< Computing standard orientation >>>"
  N = coord.shape[0]; z=[]
  # ----- <<< prepare 2D coordinate array >>>
  if d == 2:
    for i in range(N):
      z.append(coord[i,2])  
      coord[i,2] = 0.
  #----- <<< compute geometric center >>>
  center = np.mean(coord, 0)
  print "Geometric center:\n   ", center
  #----- <<< center coordinates and compute principal axes >>>
  coord = coord - center
  inertia = np.dot(coord.transpose(), coord)
  e_values, e_vectors = np.linalg.eig(inertia)
  # ----- <<< order eigenvalues (and eigenvectors) >>>
  order = np.argsort(e_values) 
  eval3, eval2, eval1 = e_values[order]
  axis3, axis2, axis1 = e_vectors[:, order].transpose()
  print "Inertia axes:\n   ", axis1, "\n   ", axis2, "\n   ", axis3
  # ----- <<< orient coordinates >>>
  Rot = np.array([axis1,axis2,axis3]) # rotation matrix
  for i in range(N):
     if d == 2: # restore z 
          coord[i,2] = z[i] 
     coord[i] = np.dot(Rot,coord[i]) # rotate
  # ----- <<< print stats >>>  
  print "XYZ axis limits:\n   " , coord[:,0].min(), "to", coord[:,0].max() 
  print "   " , coord[:,1].min(), "to", coord[:,1].max() 
  print "   " , coord[:,2].min(), "to", coord[:,2].max()  
  return(coord)

