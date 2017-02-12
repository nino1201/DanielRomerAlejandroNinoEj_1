import numpy as np
import matplotlib.pyplot as plt
L=5.
A = np.loadtxt("a.dat")
Y, X = np.mgrid[-L/2.:L/2.:512j, -L/2.:L/2.:512j]
matx=np.zeros([512,512])
maty=np.zeros([512,512])
pot=np.zeros([512,512])
for i in range(512):
    for j in range(512):
        matx[i,j]=A[i*512+j,0]
        maty[i,j]=A[i*512+j,1]
        pot[i,j]=A[i*512+j,2]
plt.figure()
plt.streamplot(X, Y, matx, maty,density=[1.5, 1.5],color=pot)    
plt.title('Campo electrico')
plt.savefig('placas.pdf')
plt.close() 
