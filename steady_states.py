import numpy as np
from scipy.integrate import odeint
from allthethings import PyPipe_ps
import matplotlib.pyplot as plt

def F(y,t,Q,L,N,D,M,a,S0,Mr):
    #(Q,L,N,D,M,a,S0,Mr) = args
    print "y = %f" %y
    p0 = PyPipe_ps(N,D,L,M,a)
    p0.Mr = Mr
    p0.S0 = S0
    S = p0.getSourceTerms(y,Q) 
    u = Q/y
    c =p0.Cgrav(y,True)
    return S/(c**2-u**2)

def main():
    L = 10
    N = 100
    M = 10
    a = 100
    D = .1
    Q =  -0.04890698 
    Mr = 0.007
    S0 = 0.042731691
    S0 = 0.116914981
    S0s = [0.042731691,0.174868713,0.116914981,0.065484615,-0.178530533,-0.017923349,0.042784691,-0.246330809]
    t0 = 0
    p0 = PyPipe_ps(N,D,L,M,a)
    y0 = p0.AofH(10,False)
    print F(y0,0,Q,L,N,D,M,a,S0,Mr)
    ts = np.linspace(0,L,100)
    '''for i in range(2):
        
        S0 = i*.05+.01
        for k in range(4):
            Mr = 0.002*k+.01
            Qt = (Q,L,N,D,M,a,S0,Mr)
            y1 = odeint(F,y0,ts,args =Qt)
            H = [p0.HofA(y,False) for y in y1]
            plt.plot(ts,H,label ="Mr = %.3f, S0=%.2f"%(Mr,S0))
            '''
    dys =[]
    for k in range(len(S0s)):
        S0 = S0s[k]
        Qt = (Q,L,N,D,M,a,S0,Mr)
        y1 = odeint(F,y0,ts,args =Qt)
        dys.append((y1[-1]-y1[0])/L)
    for k in range(len(S0s)):
        print "%f    %e"%(S0s[k],dys[k])
    #plt.plot(ts,H,label ="Mr = %.3f, S0=%.2f"%(Mr,S0))
    #plt.legend(loc='lower left')
    #plt.show()
    #r = ode(F).set_integrator('zvode', method= 'bdf', with_jacobian =False)
    #r.set_initial_value(y0,t0).set_f_params(Q)
    #dt = 0.1
    #t1 = L
    #print r.y
    #while r.successful() and r.t<t1:
    #    r.integrate(r.t+dt)
    #print("%g %g" % (r.t, r.y))

if __name__ == "__main__":
    main()


