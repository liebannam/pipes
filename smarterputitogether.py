import os
import string
import math


def getinrange(theta):   #put angles in range [-180, 180] rather than [0, 360]  (yes, in degrees)
	if theta<-180:
		theta +=360
	if theta>180:
		theta -=360
	return theta	

def colorme(x,xmin, xmax):  #return an rgb tuple for x in [xmin,xmax]. Do some kind of batshit scaling so the map is not stupid. 
	#val0 = 0 and val1 = 1
        mb = (25./255., 25./255., 112./255.) #midnight blue
        sb  = (70./255., 130./255., 180./255.) #steel blue
	iv =  (255./255., 255./255., 240./255.);#ivory
        rr =  (1,0,.0);  #red
    
        sx = (x-xmin)/(xmax-xmin)  #scaled x
	idx = 0
        delta = .1;
        cd = 1./(math.exp(-1/delta)-1);
	for i in range(len(colors)):
		if sx>=colors[i][0]:
			idx = i;
        #idx = min(idx, len(colors)-2)
	c = [1,1,1]
	for i in range(3):
		dx = sx-colors[idx][0]
		if dx>1e-5:
                        print "xmin = %f, xmax = %f x = %f idx = %d i = %d"%(xmin, xmax, x,idx,i)
			dc=colors[idx+1][i+1]-colors[idx][i+1]
		else:
			dc =0;
		#c[i] = colors[idx][i+1]+dc*dx
                w = cd*(math.exp(-sx/delta)-1); 
                c[i] = (1-w)*mb[i]+w*(iv[i]);
        print w        
        print (c[0], c[1], c[2])        	
        return (c[0],c[1],c[2])

nwrites = 199
cbarsize = 9.   #size of colorbar on actual map.  Fuck all if I know how to automatically make it look nice. I'm pretty bamboozled by all the camera angle nonsense.
cbarloc = [23,0];
tloc = [12,9];   #location of title
tscale = .8;
## values for known networks: (7deSeptiembre) cbarsize = 30, cbarloc = [10,60], tloc = [80,90],tscale = 1.5
##                            (simple3pipes)  cbarsize = 3, cbarloc = [];
rscale = 2.;  #scale the radii so that the picture has a nicer aspect ratio

fromabove = 0  #set to 0 or 1--see below
######
#cut height fields with rectangles so false color map highlights variations in h
##########
if fromabove:
    colormap = "H2Ofalse"
    colormap = "H2Oblues"
    drawpipes = 0
    shape = "cylinder"
    bckgnd = 1
#######
#cut height fields with cylinders so it looks realistic (read: amazing!) close up (only shows if pipe is full or not-pressure head not apparent)
#######
else:
    colormap = "H2Oblues"
    #colormap = "H2Ofalse"
    drawpipes = 1
    bckgnd = 1


###########
#read in network data from file created by setupandrun.cpp
###########

#mdata = open("build/output_data/mapdata.txt", 'r')
mdata = open("output_data/mapdata.txt", 'r')
nodes =[];
conns = [];
x= [];
y = [];
z = [];
linecount = 0;
for line in mdata:
        if linecount == 0:
            s = line.split();
            T = float(s[0])
	elif linecount ==1:
	    s= line.split();
	    for i in range(len(s)/2):
		    conns.append(int(s[2*i]));
		    conns.append(int(s[2*i+1]));
	elif linecount ==2:	
	    r = [rscale*float(thing) for thing in line.split()];
	    print r	
	else:
	    stuff = line.split()
	    if(len(stuff)>0):
		    nodes.append(int(stuff[0]))
		    x.append(float(stuff[1]))
		    y.append(float(stuff[2]))
		    z.append(float(stuff[3]))
        linecount +=1;
lends = conns[0::2]
rends = conns[1::2]

mdata.close()
Nnodes = len(nodes);
Nedges = len(lends);
nodetypes = [0]*Nnodes
for i in range(Nnodes):
	nodetypes[i]+=rends.count(i)
	nodetypes[i]+=lends.count(i)
print nodetypes
print "%d edges, %d nodes" %(Nedges, Nnodes)
print lends
print rends
print conns
incoming = []
for i in range(Nnodes):
	idx = [-1]
	for k in range(conns.count(i)):
		idx.append(conns.index(i,idx[k]+1))
	incoming.append([thing/2 for thing in idx[1:]])
print incoming
o_cyl = [];
i_cyl = [];
o_sph = [];
i_sph = [];
xout = [];
yout = [];
thetaout =[];
Lout = [];
hfields = [];
Ls = []
for i in range(Nedges):
	xstart = x[lends[i]]
	ystart = y[lends[i]]
	xend = x[rends[i]]
	yend = y[rends[i]]
	L = math.sqrt((xend-xstart)**2+(yend-ystart)**2);
	theta = math.degrees(math.atan2( (yend-ystart), (xend-xstart)));
	#thetas.append(theta);
	xout.append("#declare x%i = %f;\n" %(i, xstart));
	yout.append("#declare y%i = %f;\n" %( i, ystart));
	thetaout.append("#declare ang%i = %f;\n" %( i, theta));
	Lout.append("#declare L%i = %f;\n" %( i, L));
	Ls.append(L);
	print theta
	print "%f %f %f %f"% (xstart, ystart, xend, yend)


hmax = 0.
hmin = 0.

#####
#read in time step and height field maxima information 
###
sdata = open("output_data/scalings.txt", 'r')
hs = []
for line in sdata:
    stuff = line.split()
    count = int(stuff[0])
    print stuff
    hstmp = []
    for k in range(0,(len(stuff))/2):
        #hstmp.append(float(stuff[2*k]  ))
        hstmp.append(float(stuff[2*k+1]))
    print hstmp
    hmax = max(hmax, max(hstmp))
    hmin = min(hmin, min(hstmp))
    hmin = 0;
    hs.append(hstmp)
print "count =%03d hmax= %f, hmin = %f" % (count,hmax, hmin)
print hs
sdata.close()
dt = T/float(count)
print "length of hs is %d"% len(hs)

#########  
#drawing commands for the pipes
#########
for i in range(Nedges):
	ltype = nodetypes[nodes[lends[i]]]
	rtype = nodetypes[nodes[rends[i]]]
		#open at left end, closed at right
	if(ltype >1 and rtype ==1):
		o_cyl.append("cylinder{<0.,0,R%d*ro>, <1.,0,R%d*ro>, R%d*ro scale<L%d,1,1> rotate<0,0,ang%d>  translate<x%d,y%d,0> }\n" %(i,i,i,i,i,i,i))
		i_cyl.append("cylinder{<0.0,0,R%d*ri>, <.99,0,R%d*ri>, R%d*ri scale<L%d,1,1> rotate<0,0,ang%d>  translate<x%d,y%d,0> }\n" %(i,i,i,i,i,i,i))	
	#open at right end, closed at left	
	elif(ltype ==1 and rtype >1):
		o_cyl.append("cylinder{<0.,0,R%d*ro>, <1.,0,R%d*ro>, R%d*ro scale<L%d,1,1> rotate<0,0,ang%d>  translate<x%d,y%d,0> }\n" %(i,i,i,i,i,i,i))
		i_cyl.append("cylinder{<0.01,0,R%d*ri>, <1.0,0,R%d*ri>, R%d*ri scale<L%d,1,1> rotate<0,0,ang%d>  translate<x%d,y%d,0> }\n" %(i,i,i,i,i,i,i))
	#open at both ends
	else:
		o_cyl.append("cylinder{<0,0,R%d*ro>, <1.,0,R%d*ro>, R%d*ro scale<L%d,1,1> rotate<0,0,ang%d>  translate<x%d,y%d,0> }\n" %(i,i,i,i,i,i,i))
		i_cyl.append("cylinder{<0.0,0,R%d*ri>, <1.0,0,R%d*ri>, R%d*ri scale<L%d,1,1> rotate<0,0,ang%d>  translate<x%d,y%d,0> }\n" %(i,i,i,i,i,i,i))
#######
#drawing commands for the height fields
#####
	delta = r[i]/Ls[i];
	#s1 = "\n//pipe %d \nintersection{\nintersection{\nmerge{"%(i)
	if fromabove:
		s1 = "\n//pipe %d \nintersection{\nmerge{"%(i)
		shape = "box{<0,-R%d,0> <L%d,R%d,hmax%d+1> } cylinder{<0,0,0> <0,0,hmax%d+1>R%d} cylinder{<L%d,0,0><L%d,0,hmax%d+1> R%d}  }\n"%(i,i,i,i, i,i,i, i,i,i)
		sh1 = "height_field{tga fig%d rotate<90,0,0> scale<%.5f,R%d*5,hmax%d> translate<%f, R%d,-%f> scale<L%d,1,1>}}\n"%(i,1+2*delta ,i,i, -delta,i,-hmin,i)
	else:
		s1 = "\n//pipe %d \nintersection{\nintersection{\nmerge{"%(i)
		shape = "cylinder{<0,0,R%d>, <L%d,0,R%d>, R%d} sphere{<0,0,R%d> R%d} sphere{<L%d,0,R%d>R%d}  }\n"%(i,i,i, i,i,i, i,i,i)
		sh1 = "height_field{tga fig%d rotate<90,0,0> scale<%.5f,R%d*5,hmax%d> translate<%f, R%d,0> scale<L%d,1,1>}}\n"%(i,1+2*delta ,i,i, -delta,i,i)
		s1 += shape	
	#sh2 = " height_field{tga fig%d rotate<90,0,0> scale<%.5f,R%d*5,hmax%d>}\n"%(i,1 ,i,i)
	#s2 = sh1+sh2+ "translate<0,R%d,0> scale<L%d,1,1>}}\n"%(i,i)
	s2 = sh1;
###########                
#commands to cut intersections with planes so they look clean
##########
	planes =""
	for jj in range(2):
		if jj == 0:
			thisend =lends[i]
			otherend = rends[i];
			herex = x[thisend]
			sign = -1.
			where = 0.
			if(nodetypes[thisend]==2):
				close = rends.index(thisend)
				
		else:
			thisend = rends[i]
			print thisend
			otherend = lends[i]
			herex = x[thisend]
			sign = 1.
			where = Ls[i]
			if(nodetypes[thisend]==2):
				try:
					close = lends.index(thisend)
				except ValueError:
					close = rends.index(thisend)
		relevant =[ thing for thing in incoming[thisend]]; #find other pipes coming into this node
		print "i is %d and thisend is %d" %(i,thisend)
		print relevant;
		print relevant;

		if(nodetypes[thisend]==1):
			planes +="plane{ <%.1f,0,0> %.5f} \n"%(sign, where)
		elif(nodetypes[thisend] ==2):#have to figure out other pipe touching
			relevant.pop(relevant.index(i))
			neighbor = relevant[0]
			if(lends[neighbor]==thisend):
				inc = rends[neighbor];
				sign2 = -1;
			else:
				inc = lends[neighbor];
				sign2 = 1;
			theta1 = math.degrees(math.atan2(sign*(y[thisend]-y[otherend]), sign*(x[thisend]-x[otherend])))
			theta2 = math.degrees(math.atan2(sign2*(y[thisend]-y[inc]), sign2*(x[thisend]-x[inc])))
			print "thetas: %f  %f" %(theta1, theta2)
			if(sign*sign2)>0:
				angle = (theta2-theta1)-180
			else:
				angle = (theta2-theta1)
			print angle
			print getinrange(angle)
			planes+="plane{<%f,0,0>0 rotate<0,0,%f> translate <%.5f,0,0>}\n"%(sign,getinrange(angle)/2.,where) 
		elif(nodetypes[thisend] ==3):  #have to figure out other pipes touching/killmenow				
			relevant.pop(relevant.index(i))
			angles = [0,0]
			for kk in range(0,2):   #loop over neighbors
				neighbor = relevant[kk]	
				if(lends[neighbor]==thisend):
					inc = rends[neighbor];
					sign2 = -1;
				else:
					inc = lends[neighbor];
					sign2 = 1;
				print "neightbor = %i and inc = %i" %(neighbor, inc)
				print lends
				print rends
				print x
				print y
				theta1 = math.degrees(math.atan2(sign*(y[thisend]-y[otherend]), sign*(x[thisend]-x[otherend])))
				theta2 = math.degrees(math.atan2(sign2*(y[thisend]-y[inc]), sign2*(x[thisend]-x[inc])))
				print "thetas: %f  %f" %(theta1, theta2)
				if(sign*sign2)>0:
					angles[kk] = (theta2-theta1)-180
				else:
					angles[kk] = (theta2-theta1)
			print angles
			angles = [getinrange(t) for t in angles]
			print angles
			planes+="plane{<%f,0,0>0 rotate<0,0,%f> translate <%.5f,0,0>}"%(sign,angles[0]/2.,where) 
			planes+="plane{<%f,0,0>0 rotate<0,0,%f> translate <%.5f,0,0>}\n"%(sign,angles[1]/2.,where) 
#	cval = str(colorme(x,0,hmax));
	cval = 1;
	#colorme.strip(')(');
	s4 = " rotate<0,0,ang%d>  translate<x%d,y%d,0>} \n" %(i,i,i)
	hfields.append(s1+s2+planes+s4)
####
#commands to smooth out the corners/gaps with spheres at nodes
####
for i in range(Nnodes):
	if(nodetypes[i]>1):
		R = max([r[thing] for thing in incoming[1] ]);
		o_sph.append("sphere{<%.5f,%.5f,%.5f>, %.5f}\n"%(x[i], y[i], R*1.05,R*1.05));
		i_sph.append("sphere{<%.5f,%.5f,%.5f>, %.5f}\n"%(x[i], y[i], R*1.01, R*1.01));	
print "Goddamnit"


###write temporary povray plotting file with data from each time step
count = nwrites
for i in range(0,count+1):
	istring = "%03d"%i
	lines = [];
	for j in range(Nedges):
	    lines.append( "#declare fig%d = \"output_data/out%d_"%(j,j) +istring+ "\";\n")
	    lines.append("#declare R%d = %f;\n #declare hmax%d = %f;\n"%(j, r[j], j, hs[i][j]))
	fout = open("plottmp.pov", 'w')
	for j in range(0,len(lines)):
		fout.write(lines[j])
	for j in range(Nedges):
		fout.write(xout[j])
		fout.write(yout[j])
		fout.write(thetaout[j])
		fout.write(Lout[j])
	if drawpipes:
		fout.write("union{\ndifference{\nmerge{\n")
		for j in range(Nedges):
			fout.write(o_cyl[j])
		for j in range(len(o_sph)):
			fout.write(o_sph[j])
		fout.write("}")
		for j in range(Nedges):
			fout.write(i_cyl[j]);
		for j in range(len(i_sph)):
			fout.write(i_sph[j])
		fout.write(" material{ texture{pigment{ rgbf<.93,.95,.98,0.825>*0.99} finish { ambient 0.0 diffuse 0.15 reflection{0.1,0.1}specular 0.6 roughness 0.005 conserve_energy}} interior{ ior 1.6 fade_power 1001 fade_distance 0.5 fade_color <0.8,0.8,0.8> caustics 0.16}} clipped_by{plane {z 2*R0*.8}}}}\n")
	fout.write("merge{\n")
	for j in range(Nedges):
		fout.write(hfields[j])
	print 2*max(r)
        if fromabove:
            colormax = hmax-hmin;
        else:
            colormax = 2*max(r)
    	cscale = "}scale %.2f}"% (1.1*colormax);
        fout.write("texture{pigment{gradient <0,0,1> color_map {"+colormap+cscale +"finish{ambient 1.0 diffuse 0.}}}")
        
        #write the colorbar
        fout.write("#declare cbarmax = %.2f;\n" %cbarsize)
        fout.write("union {box {<0,0,0>,<1,cbarmax,0.001>pigment {gradient y color_map {"+colormap+"}scale cbarmax} finish{ambient 1.0 diffuse 0.}}")
        
        Nbars =3  #number of ticks on colorbar
        for j in range(Nbars):
            dy = float(j)*cbarsize/float(Nbars-1) 
            fout.write("box{<-.5,%.2f,0><1.5,%.2f, 0.001> pigment{rgb 1}}\n"%(dy-.1, dy+.1))
            fout.write(" text{ttf \"/usr/local/share/povray-3.7/include/timrom.ttf\" \"%.3f\" 0.15,0 pigment{rgb 1}scale %.2f translate <2.,%.2f, 0>finish{ambient 1.0 diffuse 0.}}\n"%(colormax/float(Nbars-1)*j+hmin,tscale, dy-.5))
	fout.write("translate<%.2f,%.2f,0>}\n"%(cbarloc[0], cbarloc[1]))
        fout.write(" text{ttf \"/usr/local/share/povray-3.7/include/timrom.ttf\" \"t = %.2f s\" 0.15,0 pigment{rgb 1}scale %.2f translate <%.2f,%.2f,0>finish{ambient 1.0 diffuse 0.}}\n"%(float(i)*dt,tscale, tloc[0], tloc[1]))
        fout.write("background{ rgb %d}" %bckgnd)
        fout.close()
	pngname = "tmp_"+istring+".png"
        if fromabove:
            command = "povray +A0.001 -J topview.pov Output_File_Name=movie/tmp_"+istring +".png"
        else:
            command = "povray +A0.001 -J angleview.pov Output_File_Name=movie/tmp_"+istring +".png"
        print command
	#os.system(command)

print hs
print hmax

