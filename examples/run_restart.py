import sys
sys.path.append("..")
from allthethings import PyNetwork, PyPipe_ps
from allthethings import PyBC_opt_dh
import numpy as np
import matplotlib.pyplot as plt
from writeit import rewritePipes


def idx_t(i,j,n,N):
    return (2*(N+2)*n+(N+2)*i+j)

fi = "../indata/Alamedanewer2.0.inp"
fc = "../indata/Alamedanewer2.0.config"
n0 = PyNetwork(fi,fc,1)
m_per_ft = .3048
elevs  =[417,414,324,300,275,256,232,201,177,190,192,197,199,206,209,240,252,250,283,289,316,315,334,340,334,332,333,341,387,426,417]

Ls_m = [l for l in n0.Ls]
elevs_m = [float(el)*m_per_ft for el in elevs]
Ds_m = [D/12 for D in n0.Ds]
T = 20
Mi = 1   #number of time steps in between writes
Nt = 10
Ttot = Nt*T
Np = n0.Nedges
a = 100
Ns = [int(l) for l in Ls_m]
dx = [Ls_m[i]/Ns[i] for i in range(Np)]
M = int(T*a/(max(dx)*.8))
M = (M+Mi-M%Mi)  #round it up to be an even multiple of Mi
jt = n0.nodeTypes
Nn = len(jt)
bt = [1]*Nn
bv = [0.]*Nn
r = [-1]*Nn
r[0] = 0
bt[0] =1
q0s = [0]*Np
h0s = [0.1*d for d in Ds_m]
Mrs =[0.007]*Np
print T
print M
Nstar =1 #measuring point for each pipe
Hs =np.ndarray((Np,M/Mi*Nt))

fn = "../indata/Alameda_m3"
oldinp = "../indata/Alamedanewer52.0.0.inp"
(fi, fc) = rewritePipes(fn,oldinp, Ns, Ls_m, Mrs, Ds_m, jt, bt, bv, r, h0s, q0s, T, M, a,elevs_m)
n1 = PyNetwork(fi,fc,1)
dt = n1.T/float(n1.M)
Q00 = 0.0087
for i in range(0,Np):
    d = n1.Ds[i]
    A0 = (d*d/4.)*np.ones(n1.Ns[i])
    Q0 = 0*np.ones(n1.Ns[i])
    n1.setIC(i,A0,Q0)
Qb = Q00*np.ones(M+1)
n1.setbVal(0,Qb)

Ttot= 0
Vs = [n1.getTotalVolume()]
for m in range(Nt):
    try:
        n1.runForwardProblem(dt)
    except:
        print "whoops. dt is probably too small. quitting at time T=%f"%((m-1.)*T)
        break
    for j in range(Np):
        N = n1.Ns[j]
        p0 = PyPipe_ps(N,n1.Ds[j], n1.Ls[j],M,a)
        qh = n1.qhist(j)
        Htemp = [p0.pbar(qh[idx_t(0,Nstar,n,n1.Ns[j])],False) for n in range(1,M+1,Mi)]
        Hs[j,m*(M/Mi):(M/Mi)*(m+1)] = Htemp
    Vs.append(n1.getTotalVolume()) 
    Ttot +=T
    print Ttot
    print Vs
    n1.reset()    
