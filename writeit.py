

def rewritePipes(fn,oldinp, Ns, Ls, Mrs, Ds, jt, bt, bv, r, h0s, q0s, T, M, a):
	newconfig = fn+(".config")
	newinp = fn+".inp"
	Mi = 10
#config file section titles
	Ptitle ="[PIPE_INFO]\n\
;-------------------------------------------\n\
;initial  initial\n; ID    N       h       Q\n\
;-------------------------------------------\n"
	Jtitle = "[JUNCTION_INFO]\
;-------------------------------------------------------------------------------------\n\
;{-----for junction1s-----} | {--for junction2s--} | {------for junction3s-------}|\n\
;ID   	jtype	bvaltype  bval   reflect   | offset   valveopen   | offset01   offset02   offset12 |\n\
;------------------------------------------------------------------------------------\n"
	Ttitle = "[TIME_INFO]\
;---------------------------------------\n\
;T (s)           M        Mi     a  (m/s)\n\
;-----------------------------------------\n"

#open config file and write
	with open(newconfig, 'w') as fw:

		fw.write(Ptitle)
		for j in range(len(Ns)):
			fw.write("%d    %d   %2.2f   %2.2f\n"%(j,int(Ns[j]), h0s[j],q0s[j]))
		fw.write("\n")
		fw.write(Jtitle)
		for k in range(len(jt)):
			fw.write("%d     %d     %d     %d     %d     %d     %d     %d     %d     %d\n"%( k,jt[k],bt[k],bv[k],r[k],0,1,0,0,0))
		fw.write("\n")
		fw.write(Ttitle)
		fw.write("%3.3f       %d         %d        %.1f"%(T,M, Mi, a))
	fw.close()
#now open old inp file and change pipe properties, leaving everthing else the same
	count = 0
	count2 =0

	with open(oldinp, 'r') as fold:
		with open(newinp, 'wb') as fnew:
			for line in fold:
				s = line.split()
				if '[' in line:
					fnew.write(line)
					count +=1
				if (count ==5) and (';' not in s[0]) and len(s)>1:
					fnew.write(  "\n%s               %s                	  %s                  %4.1d             %2.2f           %1.3f\n"%(s[0], s[1], s[2],Ls[count2], Ds[count2], Mrs[count2]))
					count2+=1
					if count2>=len(Ns):
						count+=1
				else:
					fnew.write(line)
	print "new files are %s and %s"%(newinp, newconfig)
	return (newinp, newconfig)

def main():
	
	fn = "test1"
	oldinp = "basecase.inp"
	Ns = [100,100,100]
	Ls = [100,100,100]
	Mrs = [0.015]*3
	Ds = [1.]*3
	jt = [1,3,1,1]
	bt = [1,1,1,1]
	bv = [0,0,0,0]
	r = [0,0,-1,-1]
	h0s = [.8,.8,.8]
	q0s = [2.,1.,1.]

	T = 18
	M = 4200
	Mi =50
	a = 100
	(fi, fc) = rewritePipes(fn,oldinp, Ns, Ls, Mrs, Ds, jt, bt, bv, r, h0s, q0s, T, M, a)
	print "fi = %s, fc = %s"%(fi,fc)

if __name__ =="__main__":	
	main()



