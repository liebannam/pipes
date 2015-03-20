'''
Use this code to
-> Run convergence test for simulated network in sim.inp and sim.config
-> Store the data in a class "conv_results", which is then pickled and read by the plotting script **plotconv.py**

Assumptions:
            -there exists a  directory called tmpconfig/ where the additional config files can be stored
            -you are happy running one of the simulations I've listed. In principle you can run any pair of .inp and .config files but I make absolutely no guarantees about plotting them in any meaningful fashion later.
            -channel.cpp and justrunit.cpp are the same as they were on February 23, 2015 (else the output formatting may be broken and what this code attempts to read from the terminal printout of channel.cpp function "quickwrite()" as called by /.justrunit will possibly fail. 
'''
import math
import subprocess
import numpy as np
import pickle


class conv_results:
    def __init__(self, X,A,H,Q,T,ymax):
        self.x = X
        self.a = A
        self.h = H
        self.q = Q
        self.t = T
        self.nref = len(X)
        self.npipes = len(X[0])
        self.hmax = ymax


def main():   
	
        npipes = 1  #number of pipes. only choose 1,2, or 3. Don't fuck with the relevant config files, please...
        nref = 5; #number of refinements, I'd recommend keeping this less than 5
        if npipes ==1:
		sim = '1pipe'   #Example from Leon (2006). Exact solution known
	if npipes ==2:
		sim = '2pipes'
	if npipes ==3:
		sim = '3pipes1' #Example illustrating triple junction sovles. Exact solution not known



#read original inp and config files
	ifile = sim+'.inp'
	cfile = sim+'.config'
	where = '../indata/'   #where original files live
	where2 = '../indata/tmpconfig/'  #where new config files live
        X = []
	A = []
	H = []
	Q = []
#create nref new config files-- double number of grid points and halve the time step from previous run
	for k in range(0,nref):
		count = 0;
		fr = open(where+cfile,'r')
		K = pow(2,k)
		newconfig = where2+sim+("%d.config"%k)
		with open(newconfig, 'w') as fw:
			for line in fr:
				#print "count = %d" %count
				#print line
				if '[' in line:
					fw.write(line)
					count +=1
				elif (count ==1) and (";" not in line) and len(line.split())>1:
					s = line.split()
					N = int(s[1]);
					fw.write(  "%s     %i	  %s     %s\n"%(s[0],K*N, s[2], s[3]))
				elif (count ==3) and (";" not in line) and len(line.split())>1:
                                        s = line.split()
					M =  int(s[1]);
					Mi = int(s[2]);
					fw.write("%s        %i        %i"%(s[0],K*M, K*Mi))
				else:
					fw.write(line)
		fr.close()
		fw.close()
		cmd = './justrunit' 
		arg = [where+ifile, newconfig]
#run the simulation as a subprocess and capture the terminal output for analysis 
		print "now running\n"+ cmd +' '+ arg[0] + ' ' +arg[1]
		print "M is %d and N is %d"%(K*M, K*N)
		proc = subprocess.Popen([cmd, arg[0] ,arg[1]], stdout=subprocess.PIPE)
		out, err = proc.communicate()
		print out
		print "C++ Runtime errors: %s \n *****" %err	
		Nn = K*N
#find relevent data within output
		x = np.zeros([npipes,Nn+1])
		a = np.zeros([npipes,Nn+1])
		h = np.zeros([npipes,Nn+1])
		q = np.zeros([npipes,Nn+1])
		ymax = 0
                wanted = out.split("x             A               h               Q") 
                print "*******\n\n"
                for j in range(npipes):
          		w = wanted[j+1].split('\n')
                    	for i in range(Nn+1):
                        	(x[j,i], a[j,i], h[j,i], q[j,i]) = [float(thing) for thing in w[i+1].split()]
				ymax = max(h[j,i],ymax)
#store all this
#eg. A = [[A0(pipe0), A0(pipe0), ..],[A1(pipe0)..]] where Ai is the ith refinement level
		X.append(x)
		A.append(a)
		H.append(h)
		Q.append(q)
        blah = wanted[0].split('\n')[-2].split(',')[0]
        T = float(blah[4:])
        ymax = round(ymax*1.05,1)
        res = conv_results(X,A,H,Q,T,ymax)
        with open('filename','wb') as fp:
            pickle.dump(res, fp)	         
        print "T = %f" %T

if __name__=="__main__":
    main()

##old and terrible way to do it
#run the simulation and store the the incoherent output to wtf.txt
# cmd = "./justrunit "+ where+inpfile +where+configfile  2>&1 | tee wtf.txt
#    os.system(cmd)

