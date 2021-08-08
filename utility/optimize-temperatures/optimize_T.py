#In this program, we estimate temperatures in a way that 
#the acceptance ratio in parallel tempering swapping becomes equal for all temperatures.
# We use the following paper for this purpose:
# Koji Hukushima, Physical Review E, 60, 3606 
# The program is written by Dr. Mojtaba Alaei.


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import sys
import argparse
parser = argparse.ArgumentParser()



print("           OPTIMIZING TEMPERATURE         ")
print("\n")


parser.add_argument('--data', '--data_file', type=str,
                    default='E.dat', help='The file contains temperature and energy')

parser.add_argument('--alpha', '--mixing_parameter', type=float,
                    default=0.5, help='Mixing parameter. If the temperatures do not converge, use a smaller value for alpha.')


parser.add_argument('--epsilon', '--convergence_threshold', type=float,
                    default=1e-6, help='Convergence threshold')

args = parser.parse_args()


dat=np.loadtxt(args.data)
T=dat[:,0]
e=dat[:,1]
E = interp1d(T, e)



def g(T1,T,T2):
    return 1/(E(T1)-E(T2)) * (1/T1*E(T1)-1/T2*E(T2)-E(T)*(1/T1-1/T2))


Tin0=dat[:,0]
N=len(Tin0)

Tout=np.zeros(N)
Tout[0]=Tin0[0]
Tout[-1]=Tin0[-1]
Tin=Tin0.copy()

eps= args.epsilon   # The convergence threshold

alpha = args.alpha

n_iter=1
while (abs(Tin-Tout).max() > eps):
   #for i in range(1,N-1):
   #   beta=1/2*((1/Tin[i])+g(Tin[i-1],Tin[i],Tin[i+1]))
   #   Tout[i]=1/beta
   # to seed up, we vectorize the three above lines in the following line:  
   Tout[1:N-1]=1/(1/2*((1/Tin[1:N-1])+g(Tin[0:N-2],Tin[1:N-1],Tin[2:N])) )
   n_iter=n_iter+1
   if (n_iter % 1000 == 0):
      print("iteration:", n_iter, "distance between Tin and Tout: %12.8f" % abs(Tin-Tout).max())
   #Tin,Tout=Tout,Tin   
   Tin[1:N-1]=(1-alpha)*Tin[1:N-1] + alpha* Tout[1:N-1]

print("\n")
print("Optimized temperatures:")
print("\n")
print("tems=", end=' ')
for T in Tout:
    print("%5.4f" % T, end=' ')

print("\n")



C0=[]
for i in range(1,N):
     C0.append((1/Tin0[i-1]-1/Tin0[i])* (E(Tin0[i-1])-E(Tin0[i])))

C=[]
for i in range(1,N):
     C.append((1/Tout[i-1]-1/Tout[i])* (E(Tout[i-1])-E(Tout[i])))

p1=plt.plot(Tin0[1:N], C0, label='initial temperatures')
p2=plt.plot(Tout[1:N], C, label='optimized temperatures')
plt.xlabel("T")
plt.ylabel(r"$\Delta E \Delta \beta$")
plt.legend()
plt.savefig("deltaE_beta.pdf")
#for i in range(1,N):
#     print("\Delta E \Delta \ beta",(1/Tout[i-1]-1/Tout[i])* (E(Tout[i-1])-E(Tout[i])))
      
