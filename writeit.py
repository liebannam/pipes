
##various utilities for dealing with .inp and .config files
#rewrite file
import numpy as np

def rewritePipes(fn, oldinp, Ns, Ls, Mrs, Ds, jt, bt, bv, r, h0s, q0s, T, M, a,elevs):
    newconfig = fn + (".config")
    newinp = fn + ".inp"
    Mi = 10
# config file section titles
    Ptitle = "[PIPE_INFO]\n\
;-------------------------------------------\n\
;initial  initial\n; ID    N       h       Q\n\
;-------------------------------------------\n"
    Jtitle = "[JUNCTION_INFO]\
;-------------------------------------------------------------------------------------\n\
;{-----for junction1s-----} | {--for junction2s--} | {------for junction3s-------}|\n\
;ID     jtype   bvaltype  bval   reflect   | offset   valveopen   | offset01   offset02   offset12 |\n\
;------------------------------------------------------------------------------------\n"
    Ttitle = "[TIME_INFO]\
;---------------------------------------\n\
;T (s)           M        Mi     a  (m/s)\n\
;-----------------------------------------\n"
    
# open config file and write
    with open(newconfig, 'w') as fw:
        fw.write(Ptitle)
        for j in range(len(Ns)):
            fw.write("%d    %d   %2.6f   %2.2f\n" %
                     (j, int(Ns[j]), h0s[j], q0s[j]))
        fw.write("\n")
        fw.write(Jtitle)
        for k in range(len(jt)):
            fw.write("%d     %d     %d     %d     %d     %d     %d     %d     %d     %d\n" % (
                k, jt[k], bt[k], bv[k], r[k], 0, 1, 0, 0, 0))
        fw.write("\n")
        fw.write(Ttitle)
        fw.write("%3.3f       %d         %d        %.1f" % (T, M, Mi, a))
    fw.close()
# now open old inp file and change pipe properties, leaving everthing else
# the same
    count = 0
    count1 = 0
    count2 = 0
    with open(oldinp, 'r') as fold:
        with open(newinp, 'wb') as fnew:
            for line in fold:
                s = line.split()
                if '[' in line:
                    fnew.write(line)
                    count += 1
                elif (count ==2) and len(s) >1 and (';' not in s[0]): 
                    fnew.write("%d %15s %2.3f %15s0 %30s ;\n"%(count1," ",elevs[count1]," "," "))
                    count1 +=1
                    if count1>=len(jt):
                        count+=1
                elif (count == 6) and len(s) > 1 and (';' not in s[0]):
                    fnew.write("%s %15s %s %15s %s %15s %4.1d %15s %2.2f %15s %1.4f\n" % \
                                (s[0]," ",s[1]," ", s[2]," ", Ls[count2]," ", Ds[count2]," ", Mrs[count2]))
                    count2 += 1
                    if count2 >= len(Ns):
                        count += 1
                else:
                    fnew.write(line)
    #print "new files are %s and %s" % (newinp, newconfig)
    return (newinp, newconfig)

def plotNetworkLayout(xs,ys,conns,ls,Np):
    from matplotlib import cm
    import matplotlib.colors as colors 
    from matplotlib import rc
    import matplotlib.pyplot as plt
    scales = [(np.sqrt(np.diff(xs[conns[k,:]])**2+ np.diff(ys[conns[k,:]])**2)/ls[k]) for k in range(Np)]
    s = np.mean(scales)
    dx  =100*s
    x1 = 10
    x2 = x1+dx
    cNorm  = colors.Normalize(vmin=0, vmax=Np+1)
    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cm.get_cmap('ocean'))
    rc('font', family='serif',size= '16')
    fig,ax = plt.subplots(figsize=(15,10))
    for k in range(Np):
        plt.plot(xs[conns[k,:]], ys[conns[k,:]], color=scalarMap.to_rgba(k), linestyle='-', linewidth=3)
    eps = x1/100.
    for j in range(len(xs)):
        x = xs[j]
        y = ys[j]
        plt.annotate( "%d "%j, xy=(x,y), xycoords='data',
        xytext=(x, y))

    plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off')
    y = 30
    plt.annotate(
        '', xy=(x1, y), xycoords='data',
        xytext=(x2,y ), textcoords='data',
        arrowprops={'arrowstyle': '|-|'})
    plt.annotate(
        "%.f m"%(dx/s), xy=(x1+1,y-1), xycoords='data',
        xytext=(x1+1, y-1), textcoords='offset points')

def getBasicConnectivity(finp):
    nodes =[]
    pipes = []
    ct = []
    xt = []
    yt = []
    ls = []
    with open(finp, 'rb') as f:
        count = 0
        for line in f:
            if '[' in line:
                p = False
            if '[COORDINATES]' in line:
                p = True
            if p:
                count+=1
                l = line.split()
                if count>2 and len(l)==3:
                    nodes.append(int(l[0]))
                    xt.append(float(l[1]))
                    yt.append(float(l[2]))
    with open(finp, 'rb') as f:
        count =0
        for line in f:
            if '[' in line:
                p = False
            if '[PIPES]' in line:
                p  =True
            if p:
                count+=1
                l = line.split()
                if count>2 and len(line)>2:
                    pipes.append(int(l[0]))
                    ct.append([int(l[1]),int(l[2])])
                    ls.append(float(l[3]))
    conns = np.array(ct)
    xs = np.array(xt)
    ys = np.array(yt)
    return (xs,ys,conns,ls)

def main():

    fn = "indata/test1"
    oldinp = "indata/3pipes3.inp"
    Ns = [100, 100, 100]
    Ls = [100, 100, 100]
    Mrs = [0.015] * 3
    Ds = [1.] * 3
    jt = [1, 3, 1, 1]
    bt = [1, 1, 1, 1]
    bv = [0, 0, 0, 0]
    r = [0, 0, -1, -1]
    h0s = [.8, .8, .8]
    q0s = [2., 1., 1.]
    elevs = [10,5,3,1]
    T = 18
    M = 4200
    Mi = 50
    a = 100
    (fi, fc) = rewritePipes(fn, oldinp, Ns, Ls,
                            Mrs, Ds, jt, bt, bv, r, h0s, q0s, T, M, a,elevs)
    #print "fi = %s, fc = %s" % (fi, fc)

if __name__ == "__main__":
    main()
