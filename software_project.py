from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def line_dir_pt(m,A):
	len = 10
	x_AB = np.zeros((3,len))
	lam_1 = np.linspace(0,10,len)
	for i in range(len):
		temp1 = A + lam_1[i]*m
		x_AB[:,i]= temp1.T
	return x_AB

#To find point of intersection of two lines
def inter(p1, m1, p2, m2):
	if abs (np.dot(p1-p2, np.cross(m1,m2))) > 1e-15 :
		return False
	s = np.dot(np.cross(p1-p2,m2),np.cross(m1,m2)) / (np.linalg.norm(np.cross(m1,m2))**2)
	return p1 - m1*s

#defining given lines : x(k) = A + k*l
A1 = np.array([1,0,0])
m1 = np.array([-1,2,2])
A2 = np.array([0,0,0])
m2 = np.array([2,-1,2]) 

def printfunc(p, m, option):
	c1 = inter(A1, m1, p, m)
	c2 = inter(A2, m2, p, m)
	y = np.linalg.norm(np.cross(m, np.cross(m1, m2))) == 0
	if isinstance (c1,bool) or isinstance (c2,bool) or y == False:
		print("option ", option, " does not satisfy the condition") 
	else:
		print("option ", option, " satisfies the given condition")

#setting up plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

#Lines to be verified 
A3_A = (2/9)*np.array([4,1,1])
m3_A = np.array([2,2,-1])
A3_B = (2/9)*np.array([2,-1,2])
m3_B =  np.array([2,2,-1])
A3_C = (1/3)*np.array([2,0,1])
m3_C = np.array([2,2,-1])
A3_D = np.array([0,0,0])
m3_D = np.array([2,2,-1])

#generating points in line 
l1_p = line_dir_pt(m1,A1)
l2_p = line_dir_pt(m2,A2)
l3_A = line_dir_pt(m3_A,A3_A)
l3_B = line_dir_pt(m3_B,A3_B)
l3_C = line_dir_pt(m3_C,A3_C)
l3_D = line_dir_pt(m3_D,A3_D)

printfunc(A3_A, m3_A,"A")
printfunc(A3_B, m3_B,"B")
printfunc(A3_C, m3_C,"C")
printfunc(A3_D, m3_D,"D")

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
