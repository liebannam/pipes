
import time
import numpy as np
import pickle
import sys
sys.path.append('/Users/anna/anaconda/lib/python2.7/site-packages/')
import networkx as nx
from scipy import stats
#from networkx import graphviz_layout

from matplotlib import rc
from matplotlib.mlab import find
from  matplotlib import pyplot as plt
rc('font', family='serif',size= '16')

sys.path.append("..")
from allthethings import PyNetwork, PyPipe_ps
from allthethings import PyBC_opt_dh
from writeit import writePipes
import os



def getCoordinates(conns, scale=100):
    Np=np.shape(conns)[0]
    xs = []
    ys = []
    G=nx.Graph()
    for k in range(Np):
        G.add_edge('%d'%conns[k,0],'%d'%conns[k,1])
    #pos=nx.graphviz_layout(G) # positions for all nodes
    pos = nx.spring_layout(G)
    for k in range(len(G.nodes())):
        s = '%d'%k
        xs.append(scale*pos[s][0])
        ys.append(scale*pos[s][1])
    return (xs,ys)

def getNodeTypes(conns):
    nl = conns[:,0]
    nr = conns[:,1]
    return [len(find(nl==k))+len(find(nr==k)) for k in range(len(np.unique(conns)))]

def drawNetwork(conns, place):
    G=nx.Graph()
    Np=np.shape(conns)[0]
    for k in range(Np):
        G.add_edge('%d'%conns[k,0],'%d'%conns[k,1])
    pos=nx.graphviz_layout(G) # positions for all nodes
        # nodes 
    plt.clf()
    plt.figure(1)
    nx.draw_networkx_nodes(G,pos,node_size=300,node_color = '#A0CBE2', iter=500)
    nx.draw_networkx_edges(G,pos,width=4,edge_cmap=plt.cm.coolwarm)
    nx.draw_networkx_edges(G,pos,width=4.0,alpha=0.3)
        # labels
    nx.draw_networkx_labels(G,pos,font_size=8,font_family='sans-serif')
    plt.axis('off')
    plt.savefig(place+'/layout.png')


def writeFiles(case,fn):
    if case==0:
        '''this is from Leon (2006) (pg 805)'''
        conns = np.array([[0,1]])
        nodeTypes = getNodeTypes(conns)
        Np = np.shape(conns)[0]
        Nn = len(np.unique(conns))
        Ls = [1000]
        Ns = [100]
        Mrs = [0.0]
        Ds = [2.5]
        jt = nodeTypes
        #reflect everything
        bt = [1,1]
        r = [0,1]
        bv = [2,0]
        h0s = [0.5]
        q0s = [2]
        T = 300
        a = 10
        elevs = [0]*Nn
        descrip = "pipe-filling bore in sewer flow from Leon (2006) pg 805"
    if case==1:
        conns = np.array([[0,1],[1,2]])
        nodeTypes = getNodeTypes(conns)
        Np = shape(conns)[0]
        Nn = len(unique(conns))
        Ls = [10,10]
        Ns = [Ls[k]*4 for k in range(Np)]
        Mrs = [0.0]*Np
        Ds = [0.1]*Np
        jt = nodeTypes
        bt = [1,1,1]
        #reflect everything 
        r = [1,1,1]
        bv = [0,0]
        h0s = [6,.01]
        q0s = [0.01,0]
        T = 30
        a = 100
        elevs = [0,0,0]# this gives weird crappy peak in spatial data
        descrip = "dam-break problem in two pipes--unstable for linear evaluation of A*"
    
    elif case ==3:
        conns = np.array([[0,1],[1,2],[1,3],[2,4], [3,4],[4,5]])
        nodeTypes = getNodeTypes(conns)
        Np = np.shape(conns)[0]
        Nn = len(np.unique(conns))
        Ls = [50,150,25,25,25,50]
        Ns = [Ls[k]*2 for k in range(Np)]
        Mrs = [0.007]*Np
        Ds = [0.1]*Np
        jt = nodeTypes
        bt = [1]*Nn
        #reflect everything except at inlet node
        r = [1]*Nn
        r[0]=0
        bv = [0]*Nn
        #specify pressure here
        bv[0] = 0.015473914995245055
        bt[0]=0
        h0s = [0.01]*Np
        h0s[0]=1
        #h0s = [1.,2.,3.,4.,5.,1.]
        q0s = [0.00]*Np
        T = 600
        a = 10
        elevs = [0]*Nn
        descrip= "pressurizing network with loop. specify pressure at inlet (nonconstant)"
    
    (xs,ys) = getCoordinates(conns)
    dx = [Ls[i]/float(Ns[i]) for i in range(Np)] 
    M = int(T*a/(max(dx)*.8))
    (fi, fc) = writePipes(fn,conns, xs,ys, Ns, Ls, Mrs, Ds, jt, bt, bv, r, h0s, q0s, T, M, a, elevs)
    return (fi,fc,descrip,conns)

def writeSummaryInfo(case,place,n0,V0,Vf, simtime, descrip):
    if case==0:
        N = n0.Ns[0]
        A =n0.q(0)[0:N]
        Q = n0.q(0)[N:]
        T = n0.T
        where = find(abs(np.diff(A))>1e-1)
        lower =where[0]
        upper = where[-1]
        L = n0.Ls[0]
        dx = L/float(N)
        s1 = (lower*dx-L)/T
        s2 = (upper*dx-L)/T  
        Al = A[0]
        Ar = A[-1]
        Ql = Q[0]
        Qr = Q[-1]

        loc = n0.Ls[0]+((Ql-Qr)/(Al-Ar))*n0.T
        x = np.linspace(0,L,N)
        y = np.linspace(0,1.2)
        plt.plot(x,n0.getHofA(0),loc*np.ones(np.size(y)),y,'k:')
        ax = plt.gca()
        ax.set_ylim(0,1.2)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('water height (m)')
        ax.set_title('Water profile at time T= %.f s'%T)
        print "shock speed based on Rankine Hugoniot is %.3f m/s"%((Ql-Qr)/(Al-Ar))
        print "numerical shock speed is %.3f +/- %.3f m/s"%((s1+s2)/2, abs(s1-s2))
        print "unnacounted for volume V(T)-(Q(0,T)*T+V(0))= %e m^3"%(Vf-V0-2*T)
        print "compute time is %f s"% (simtime)
        plt.savefig(place+'/results.eps')
        with open(place+"/info.txt",'wb') as f:
            f.write("Case %d\n"%case)
            f.write(descrip+'\n')
            f.write("Code is at state...\n")
            f.write("shock speed based on Rankine Hugoniot is %.3f m/s\n"%((Ql-Qr)/(Al-Ar)))
            f.write("numerical shock speed is %.3f +/- %.3f m/s\n"%((s1+s2)/2, abs(s1-s2)))
            f.write("unnacounted for volume V(T)-(Q(0,T)*T+V(0))= %e m^3\n"%(Vf-V0-2*T))
            f.write("compute time is %f s\n"% (simtime))
            f.write("At time T=%f s\n cell A  Q\n"%T)
            for k in range(N):
                f.write("%d   %f   %f\n"%(k,A[k], Q[k]))
    elif case ==1:
        import matplotlib.colors as colors  
        KK = 100
        levels = range(0,KK+1,1)
        CS3 = plt.contourf( [[0,0],[0,0]], levels, cmap=cm.get_cmap('ocean'))#fake contour plot to get colorbar
        plt.clf()
        cNorm  = colors.Normalize(vmin=0, vmax=KK+1)
        fig, ax = plt.subplots(nrows=2, figsize = (16,8))
        x = [np.linspace(0,Ls[0],Ns[0]),np.linspace(Ls[0],2*Ls[0],Ns[0])]
        scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cm.get_cmap('ocean'))
        for i in range(2):
            for k in range(0,KK):
                ax[0].plot(x[i],n0.pressureSpaceSeries(i,k), c= scalarMap.to_rgba(k), label='t=%.2f'%(dt*10*k))
                ax[1].plot(x[i],n0.pressureSpaceSeries(i,M-80*k), c= scalarMap.to_rgba(k), label='t=%.2f'%(dt*10*k))       
                ax[i].set_ylabel('P (m)')

        ax[0].set_title('%.3f to %.3f seconds'%(0,k*dt))
        ax[1].set_title('%.3f to %.3f seconds'%((M-KK+1)*dt,M*dt))
        ax[1].set_ylim(0,.05)    
        ax[0].set_xticklabels([])
        ax[1].set_xlabel('x (m)')
        cbaxes = fig.add_axes([.91, 0.1, 0.03, 0.8])
        plt.suptitle('pressure as a function of x, 100 snapshots in time')
        cb=fig.colorbar(CS3,cax=cbaxes)
        plt.savefig(place+'/results.eps')
        print "unnacounted for volume V(T)-(Q(0,T)*T+V(0))= %e m^3"%(Vf-V0)
        print "compute time is %f s"% (simtime)
        N = n0.Ns[0]
        A = np.concatenate(n0.q(0)[0:N],n0.q(1)[0:N])
        Q = np.concatenate(n0.q(0)[N:2*N],n0.q(1)[N:2*N])
        with open(place+"/info.txt",'wb') as f:
            f.write("Case %d\n"%case)
            f.write(descrip+'\n')
            f.write("Code is at state...\n")
            f.write("unnacounted for volume V(T)-V(0))= %e m^3\n"%(Vf-V0))
            f.write("compute time is %f s\n"% (simtime))
            f.write("At time T=%f s\n cell A  Q\n"%T)
            for k in range(N):
                f.write("%d   %f   %f\n"%(k,A[k], Q[k]))

    elif case==3:
        N = n0.Ns[0]
        dt = n0.T/float(n0.M)
        A =n0.q(0)[0:N]
        Q = n0.q(0)[N:]
        t = np.linspace(0,n0.T,n0.M+1)
        Pmax = 0
        Np = len(n0.Ns)
        fig, ax = plt.subplots(nrows=2, figsize = (10,10))
        xx = np.linspace(0,n0.T,n0.M+1)
        Q0 = (1+0.2*np.exp(-(xx-n0.T/2)**2/(n0.T/10)**2))*0.008

        print "pipe  T_p(s)"
        for i in range(0,Np):
            P = n0.pressureTimeSeries(i,n0.Ns[i]/2)
            print"%d     %.1f" %(i,(find(P>.1)[0])*dt)
            Pmax = max(max(P),Pmax)
            ax[0].plot(t,P, label='%d'%i)
        ax[0].set_title('pressure time series in middle of each pipe')
        ax[0].legend(loc = 'upper right')
        ax[1].plot(xx,Q0)
        ax[1].set_yticks([0,0.01])
        ax[0].set_xticklabels([])
        ax[1].set_ylabel('boundary value for A (m^2)')
        ax[0].set_ylabel('pressure (m)')
        ax[1].set_xlabel('time (s)')
        print "max pressure is %f m"%Pmax
        plt.savefig(place+'/results.eps')
        with open(place+"/info.txt",'wb') as f:
            f.write("Case %d\n"%case)
            f.write(descrip+'\n')
            f.write("Code is at state...\n")
            f.write("pressurization time T_p for middle of each pipe\n")
            f.write("pipe  T_p(s)\n")
            for i in range(0,Np):
                 P = n0.pressureTimeSeries(i,n0.Ns[i]/2)
                 f.write("%d     %.1f\n" %(i,(find(P>.1)[0])*dt))
            f.write("max pressure is %f m\n"%Pmax)


def runSim(case,fn,place,descrip):
    (fi,fc,descrip,conns) = writeFiles(case,fn)
    t0= time.clock()
    n0=PyNetwork(fi,fc,1)
    dt = n0.T/(float(n0.M))
    V0 = n0.getTotalVolume()
    if case==0:
        n0.runForwardProblem(dt)
    elif case==1:
        n0.runForwardProblem(dt)
    elif case ==3:
        T = n0.T
        xx = np.linspace(0,T,n0.M+1)
        Q0 = (1+0.2*np.exp(-(xx-T/2)**2/(T/10)**2))*0.008
        n0.setbVal(0,Q0)
        n0.runForwardProblem(dt)
    Vf = n0.getTotalVolume()
    tf = time.clock()
    simtime = tf-t0
    writeSummaryInfo(case,place,n0,V0,Vf,simtime, descrip)

def main():
    #available cases: 0 (one pipe) 1 (2 pipes) and 3 (loop with 6 pipes).
    #pending cases: 1 (three pipes, diverging), 2 (three pipes, converging), 4 (branching tree) 
    case = 1
    fn = '../indata/sanity_check_case%d'%case
    #make place to store data if it doesn't exist yet
    if not os.path.exists("sanity_check/"):
        os.makedirs("sanity_check")
    place = r'sanity_check/case%d'%case
    if not os.path.exists(place):
        os.makedirs(place)
    #write files and run simulation
    (fi,fc,descrip,conns) = writeFiles(case,fn)
    print "\n%s \n"%('*'*20)
    print "Case %d"%case
    runSim(case,fn,place, descrip)
    print "Case %d"%case
    print "\n%s \n"%('*'*20)
    print descrip
    drawNetwork(conns, place)
    print "further info written to %s"%(place+"/info.txt")
    print "sim results plotted in %s" %(place+'/results.eps')
    print "network layout plotted in %s"%(place+'/layout.png')
    print "\n%s \n"%('*'*20)


if __name__=='__main__':
    main()



