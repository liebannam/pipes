import math
from scipy import integrate, optimize
import numpy as np
import  matplotlib.pyplot as plt
import numpy. polynomial.chebyshev as cheb
from matplotlib import rc
from  interp_phi import *

def main():
	pi = np.pi;
	N = 20
	a1 = pi/8.;  
	a2 = pi/4.;
	p1 = 1./3.
	p2 = 5./12.
    	x = getChebNodes(N)
	ax1 = fA(x,p1,0.,a1);
	ax2 = fA(x,p2,a2, a1);
	theta1 = getTheta(ax1)
	theta2 = getTheta(ax2)
	phi1= phi(theta1)
	phi2 = phi(theta2)
	fp1 = cheb.chebfit(x,phi1,N)
	fp2 = cheb.chebfit(x,phi2,N)
#	print fp1
#	print fp2
	print "Disagreement at pi/8 is %.15f" %(cheb.chebval(1,fp1)-cheb.chebval(1,fp2));
	phim1 = cheb.chebval(1,fp2);
	phim2 = cheb.chebval(-1,fp2);
#	print phim1;
#	print phim2;

	K = 4;
	M = 3*pow(2,K);
	M = 60;
	K = 5;
	alphas = np.array([i/float(M) for i in range(2*K,M+3)]);
	Atrue1 = np.linspace(0,pi/8.,200)
	Atrue2 = np.linspace(pi/8.,pi/4.,200)
	phis1 = cheb.chebval(fX(Atrue1, p1,0,a1),fp1)
	phis2 = cheb.chebval(fX(Atrue2, p2,a2,a1),fp2)
	#alphas = [1./3.]
#	print cheb.chebval(fX(0,p1,0,a1),fp1)
#	print cheb.chebval(fX(pi/8,p1,0,a1),fp1)-phim1
	r = np.zeros((len(alphas),2))
	err = np.zeros((len(alphas),2))
	for k in range(len(alphas)):
		p= alphas[k];
		sx1 = fA(x,p,0,phim1);
        	sx2 = fA(x,p,phim2, phim1);
		A1 =np.zeros(len(sx1));
		A2 =np.zeros(len(sx1));
	#	print "p = %f" %p
	#	print sx1
	#	print sx2
		for i in range(1,len(sx1)-1):
			A1[i] = optimize.ridder(lambda x: -cheb.chebval(fX(x,p1,0,a1),fp1)+sx1[i],0,pi/8.+.1);
		A1[-1] = pi/8.;
		A2 = [optimize.ridder(lambda x: cheb.chebval(fX(x,p2,a2,a1),fp2)-sxi,pi/8., pi/4.) for sxi in sx2];
	#	print A1
	#	print A2
		fa1 = cheb.chebfit(x,A1,N)
  		fa2 = cheb.chebfit(x,A2,N)
	#	print fa1
	#	print fa2
		r[k,0] = abs(fa1[-1])
        	r[k,1] = abs(fa2[-1])
        	err[k,0] =np.linalg.norm(Atrue1-cheb.chebval(fX(phis1,p,0,phim1),fa1))
        	err[k,1] =np.linalg.norm(Atrue2-cheb.chebval(fX(phis2,p,phim2,phim1),fa2))
	print "alpha      r1      r2     err1      err2"
	for j in range(len(r)):
		print " %f   %e   %e   %e   %e" %(alphas[j],r[j,0], r[j,1],err[j,0], err[j,1])   
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')

	plt.subplots_adjust(hspace=0.4)
	plt.subplot(211)
	plt.semilogy(alphas,r)
	plt.legend([r"$\phi^{-1}_1$", r"$\phi^{-1}_2$"])
	xlab =  ['1/6','1/4', '1/3', '2/5','1/2','3/5', '2/3','3/4', '5/6', '1']
	ticks = [1./6., 1./4., 1./3.,2./5., 1./2.,3./5., 2./3.,3./4., 5./6., 1.];
	plt.xticks(ticks, xlab)
	plt.grid(True)
	plt.title(r'$\phi^{-1}(x) \approx \sum_k a_k T_k(cx^{\alpha}-1)$')
	#plt.ylabel(r'a_N')
	plt.ylabel('coefficient decay')
	plt.xlabel(r'$\alpha$')
	plt.subplot(212)
	plt.xlabel(r'$\alpha$')
	plt.ylabel('Error')
	#plt.ylabel(r'$\phi^{-1}_i(x)-\tilde{\phi}_i^{-1}(x)$')
	plt.semilogy(alphas, err)
	plt.xticks(ticks, xlab)
	plt.grid(True)
	plt.title('Error')
	plt.legend([r"$\phi^{-1}_1$", r"$\phi^{-1}_2$"])
	#plt.show()
	p3 = 1.
	p4 = 3./5;
	sx1 = fA(x,p3,0,phim1);
	sx2 = fA(x,p4,phim2, phim1);
	A1 =np.zeros(len(sx1));
	A2 =np.zeros(len(sx1));
#	print sx1
#	print sx2
	for i in range(1,len(sx1)-1):
		A1[i] = optimize.ridder(lambda x: -cheb.chebval(fX(x,p1,0,a1),fp1)+sx1[i],0,pi/8.+.1);
	A1[-1] = pi/8.;
	A2 = [optimize.ridder(lambda x: cheb.chebval(fX(x,p2,a2,a1),fp2)-sxi,pi/8., pi/4.) for sxi in sx2];
#	print A1
#	print A2
	fa1 = cheb.chebfit(x,A1,N)
	fa2 = cheb.chebfit(x,A2,N)
#	print fa1
#	print fa2
	
	print "Power for phi1 is %f"% p1
	print "Power for phi2 is %f"% p2
	print "Power for phi1^(-1) is %f"% p3
	print "Power for phi2^(-1) is %f"% p4
	fnames = ["phiofA1.txt", "phiofA2.txt","Aofphi1.txt", "Aofphi2.txt"]
	coeffs = np.zeros((N+1,4))
	coeffs[:,0] = fp1;
	coeffs[:,1] = fp2;
	coeffs[:,2] = fa1;
	coeffs[:,3] = fa2;
	for i in range(4):
		f = open(fnames[i],"w")
	        for j in range(len(fa1)):
			f.write("%.16f, " %(coeffs[j,i]) ) 
	    	f.close()
	for i in range(N+1):
		print "%.16f,  " %x[i]


if __name__ == "__main__":
	main()


