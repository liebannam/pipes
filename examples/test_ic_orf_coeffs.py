import sys
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib import rc
rc('font', family='serif',size= '16')
sys.path.append("..")
from allthethings import PyNetwork, PyPipe_ps
from allthethings import PyBC_opt_dh
from writeit import *

def idx_t(i,j,n,N):
    return (2*(N+2)*n+(N+2)*i+j)


Ls = [700,300,100,1000]
conn =np.array([[0,1],[1,2],[1,3],[2,4]])
Np = len(Ls)
T = 100  #set simulation time per run
Mi = 10   #number of time steps in between writes
Nt =20    #number of simulations to do
Ttot = Nt*T  # total simulation time
Np = len(Ls)#n0.Nedges
a = 50    #pressure wave speed (this determines slot width)

Ns = [int(l)/2 for l in Ls]   # set dx = 1 meter by assigning 1 grid cell per meter of length
Ds = [0.1]*Np
dx = [Ls[i]/Ns[i] for i in range(Np)]  
M = int(T*a/(max(dx)*.8))*4#set time steps based on dx to assure CFL condition
M = (M+Mi-M%Mi)  #round it up to be an even multiple of Mi

dt = T/float(M)

jt = [1,3,1,2,1]
Nn = len(jt)
r =[0,1,0,0,0]
bt = [1,1,2,1,2]
bv = [0,0,0,0,0]
h0s = [0,0,0,0]
q0s = [0,0,0,0]
elevs = [130,60,60,60,125]
Mrs = [0.015]*Np
xs = [200,100,70,75,0]
ys = [0,100,130,150,180]
plt.plot(xs,ys)
#create matrices to store pressure and velocity data
# i^th row of Hs is the i^th timeslice, where we take Mi time steps of length dt in between storing each time slice
# the i^th time slice consists of the [H0[i], H1[i], ...H_Np[i]], 
# where H0[i] is the pressure data as a function of x in pipe 0 at time slice i, etc
# same indexing scheme holds for Us
Hs =np.zeros((M/Mi*Nt,sum(Ns)))
Us =np.zeros((M/Mi*Nt,sum(Ns)))#

#new prefix for .inp and .config files that will be used for actual runtime parameters
fn = "../indata/Alameda_short"

#write new files [fn].inp and [fn].config with new info 
#(fi, fc) = rewritePipes(fn,oldinp, Ns, Ls, Mrs, Ds, jt, bt, bv, r, h0s, q0s, T, M, a, elevs)
(fi, fc) = writePipes(fn,conn, xs,ys, Ns, Ls, Mrs, Ds, jt, bt, bv, r, h0s, q0s, T, M, a, elevs)
m2psi = 1.42 #conversion factor, meters to psi
m32gal=264.172052 #conversion factor m^3 to gallons
for J in range(10):
    for I in range(10):
        t0 = time.clock()
        h0 = 0.05*J*Ds[0]+0.005
        orf = 0.05*I*Ds[0]+0.005
        n1 = PyNetwork(fi,fc,1)
        p0 = PyPipe_ps(Ns[0],Ds[0],Ls[0],M,a)
        A0 = p0.AofH(h0, True)
        Q0 = 0
        Qin = 0.0087
        Qout = 0.001
        Qin_t=Qin*np.ones(M+1)
        z = orf
        o = np.ones(Ns[0])
        n1.setIC(1,A0*o,Q0*o)
        n1.setbVal(0,Qin_t)
        n1.setbVal(1,z*np.ones(M+1))
        n1.setbVal(2,z*np.ones(M+1))


        t0= time.clock() 
        Vs = [n1.getTotalVolume()]  #store the system volume at times [0, T, 2*T, ...Nt*T]
        for m in range(Nt):
            n1.runForwardProblem(dt)#run the forward problem
            print 'simulation time =%f s'%(T*(m+1))#show what time we're at
            Ntot = 0
            for j in range(Np):
                N = n1.Ns[j]
                qh = n1.qhist(j)
                for n in range(1,M+1,Mi):
                    Px = n1.pressureSpaceSeries(j,n)
                    Utemp = [qh[idx_t(1,k,n,N)]/qh[idx_t(0,k,n,N)] for k in range(1,N+1)]
                    Hs[(n-1)/Mi+m*(M/Mi),Ntot:Ntot+N] = Px
                    Us[(n-1)/Mi+m*(M/Mi),Ntot:Ntot+N] = Utemp
                Ntot+=N
            Vs.append(n1.getTotalVolume()) 
            print m
            print Vs
            n1.reset()  #reset internal counter to zero so network can be run again
        elapsed_t = time.clock()-t0
        print "Wall clock time = %f"%elapsed_t
        df = pd.DataFrame(data=Hs.astype(float))
        df.to_csv('../../../Desktop/filling_sims/data_orf_%d_h0_%d_I%dJ%d.csv'%(int(np.ceil(orf/Ds[0]*100)),h0,I,J), float_format='%.4f')

        fig,ax = plt.subplots(figsize = (15,5))
        what = [25]*Np
        interesting = np.arange(0,Np)
        t = np.linspace(0,Ttot, M/Mi*Nt)
        cNorm  = colors.Normalize(vmin=0, vmax=Np+1)
        smap = cm.ScalarMappable(norm=cNorm, cmap=cm.get_cmap('ocean') )
        for k in interesting:
            ax.plot(t/60,m2psi*Hs[:,sum(Ns[0:k])+what[k]],label="pipe %d"%k,color = smap.to_rgba(k),lw=2)    
        ax.set_xlabel('t (min)')
        ax.set_ylabel('pressure head (psi)')
        plt.legend(loc = 'upper left')
        plt.gca().xaxis.grid(True)
        plt.savefig('../../../Desktop/filling_sims/data_orf_%d_h0_%d_I%dJ%d.eps'%(int(np.ceil(orf/Ds[0]*100)),h0,I,J))
        print "h0 = %f and orf = %f and inflow volume = %f"%(h0, orf, Vs[-1]-Vs[0])
