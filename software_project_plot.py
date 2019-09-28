from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

#setting up plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

l1_p = np.loadtxt('l1_p.dat',dtype='double')
l2_p = np.loadtxt('l2_p.dat',dtype='double')
l3_A = np.loadtxt('l3_A.dat',dtype='double')
l3_B = np.loadtxt('l3_B.dat',dtype='double')
l3_C = np.loadtxt('l3_C.dat',dtype='double')
l3_D = np.loadtxt('l3_D.dat',dtype='double')

#plotting line
plt.plot(l1_p[0,:],l1_p[1,:],l1_p[2,:],label="Line L1")
plt.plot(l2_p[0,:],l2_p[1,:],l2_p[2,:],label="Line L2")
plt.plot(l3_A[0,:],l3_A[1,:],l3_A[2,:],label="Option A")
plt.plot(l3_B[0,:],l3_B[1,:],l3_B[2,:],label="Option B")
plt.plot(l3_C[0,:],l3_C[1,:],l3_C[2,:],label="Option C")
plt.plot(l3_D[0,:],l3_D[1,:],l3_D[2,:],label="Option D")

#show plot
plt.xlabel('$x$');plt.ylabel('$y$')
plt.legend(loc='best');plt.grid()
plt.show()
