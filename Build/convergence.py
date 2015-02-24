'''
Run convergence test for simulated network in sim.inp and sim.config
Assumptions: -.inp and .config names match
             - data is stored in outdata/
             -
'''
import math
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import pylab as py

def main():
   
#read original inp and config files
	npipes = 1

	if npipes ==1:
		sim = '1pipe'   #Example from Leon (2006). Exact solution known
	if npipes ==3:
		sim = '3pipes1' #Example illustrating triple junction sovles. Exact solution not know
	
	ifile = sim+'.inp'
	cfile = sim+'.config'
	where = '../indata/'   #where original files live
	where2 = '../indata/tmpconfig/'  #where new config files live
	nref = 5; #number of refinements
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
                
                print len(wanted)
                print wanted[1]
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
	ymax = round(ymax*1.05,1)
#plot the results
	cNorm = colors.Normalize(vmin = 0, vmax = nref+1)
        scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = 'YlOrRd')
        fig, ax  = plt.subplots(nrows=max(npipes,2))
        legs = []
        if npipes>1:
		explain = ['to pipes 1 and 2', 'outflow', 'reflection']
        else:
		explain = ['reflection']
        for i in range(nref):
		colorVal = scalarMap.to_rgba(nref-i)
            	dx = X[i][0][1]-X[i][0][0];
            	legs.append("dx=%2.2f"%dx)
            	for jj in range(npipes):
			j = npipes-jj-1
                	ax[j].plot(X[i][j], H[i][j], color = colorVal)
                	ax[j].set_ylim([0,ymax])
                	ax[j].set_axis_bgcolor('black')
                	ax[j].set_yticks(np.arange(0,5)/(4.)*ymax)
                	ax[j].set_ylabel('pipe %d'%j)
			ax[j].tick_params(axis='y', colors='.5')
			ax[j].tick_params(axis='x', colors='.5')
                	ax[j].set_xticklabels([])
 		if npipes ==1:
			Ni = N*pow(2,i)+1
			htrue = np.zeros(Ni)
			for ii in range(Ni):
				if ii/(float(Ni))<.5828:  #using R-H conditions to compute shock speed
					htrue[ii] =  .5
				else:
					htrue[ii] =1.122893 
			ax[1].plot(X[i][0], H[i][0]-htrue, color = colorVal)
			ax[1].set_ylabel('error')
			ax[1].set_axis_bgcolor('black')
			ax[1].tick_params(axis='y', colors='.5')
			ax[1].tick_params(axis='x', colors='.5')
			ax[0].set_ylabel('h')
#failed legend attempts...bloody hell, Matplotlib...
	#plt.legend(legs,bbox_to_anchor=(0, 1.02, 1, .1), ncol = 5)
        #plt.legend(legs,bbox_to_anchor=(0., 1.02, 1,.102), loc=9, ncol=nref, mode="expand", borderaxespad=0.)
        #plt.legend(legs,bbox_to_anchor=(0.5,1.1), loc=9, ncol=nref, mode="expand", borderaxespad=0.)
	#plt.legend(legs,bbox_to_anchor=(.5, 3.5),loc ='center',ncol = nref,  mode="expand",borderaxespad=0.)
        #plt.legend(legs, 'upper left')
        
	box = ax[-1].get_position()
	ax[-1].set_position([box.x0, box.y0 + box.height * 0.1,box.width, box.height * 0.9])
# Put a legend below current axis
	ax[-1].legend(legs,loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=5)
	blah = wanted[0].split('\n')[-2].split(',')[0]
        fig.suptitle('Water height along pipes at time %s '%blah, fontsize =10)
        leg = plt.gca().get_legend()
        legt = leg.get_texts()
        frame = leg.get_frame()
        plt.setp(legt, fontsize=10, color = 'white')
        frame.set_facecolor('black') 
	ax[npipes-1].set_xlabel('x (m)')
        for j in range(npipes):
		ax2 = ax[j].twinx()
            	ax2.set_yticklabels([])
            	ax2.set_xticklabels([])
            	ax2.set_ylabel(explain[j], fontsize = 10)
	    	ax2.tick_params(axis='y', colors='white')
        #print out
	figname ='conv_results.png' 
	py.savefig(figname)
	print "figure saved as %s" %figname
	        

#give up and work at cheeseboard 

if __name__=="__main__":
    main()

##old and terrible way to do it
#run the simulation and store the the incoherent output to wtf.txt
# cmd = "./justrunit "+ where+inpfile +where+configfile  2>&1 | tee wtf.txt
#    os.system(cmd)

