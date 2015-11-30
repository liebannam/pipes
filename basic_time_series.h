/*This file is part of Pipes.

    Pipes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Pipes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Pipes.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef BASIC_TIME_SERIES_H
#define BASIC_TIME_SERIES_H


void getTimeSeries(vector<Real> & bvals, vector<Real> &x, const int m, const int M, double T, int Fourier)
       	
//bvals is time series of M+1 values at t =0, T/M, ...T
//x is a length m vector. either contains Fourier modes (Fourier==1) or of Hermite spline components (Fourier ==0)
//we have m/2 values in Fourier series or m/2 spline interpolation points. Interpolation dt = 1/(m/2-1)
//
{
	if(Fourier)	//Discrete Fourier modes:
	//x[k] = a_k (k<=m/2)
	//x[k] = b_j (where j =k-m/2, for k>m/2)
	//      x =     [a_0, a_1, ...a_m/2, b_1,   b2   ...b_(m/2-1)] 
	//                |     |       |     |     |
	//index, k        0     1      m/2  m/2+1  m/2+2       m-1     
	{
		double t;
		for(int nn = 0; nn<M+1; nn++)
		{
			bvals[nn] = 0.5*x[0];
			t = (double)(nn)/(double)M;
			for (int k = 1; k<m/2; k++)
			{
				bvals[nn] += x[k]*cos(2.*PI*(double)k*t) + x[k+m/2]*sin(2.*PI*(double)k*t);
		//		cout<<k<< " "<< x[k+m/2+1]<<endl;
			}
			bvals[nn] +=0.5*x[m/2]*cos(PI*(double)m*t);
		//	printf("n = %d, t/T = %f, gah %f\n", n, t, bvals[n]);
		//	cout<<"!!!!!"<<bvals[n]<<endl;
		}
	}
	else//Hermite Spline
	//Hermite spline evaluation, adapted from C. Rycroft's code. 
	// Evaluates spline x at point t, where
	// x = [f(t0), f'(t0), f(t1), f'(t1), ...f(t_m), f'(tm)]
	//      |       |        |      |          |       |
	// i=   0       1        2      3          2*m    2*m+1
	// x0 = 0, xm = T

	{
		double t, Dt,dt;
		dt = T/M;//time series dt
		Dt = T/((double)m/2.-1.);//spline dt. Careful! This is NOT the dt in your problem!
		for(int nn = 0; nn<M+1; nn++)
		{
			t= nn*dt;
			int j = int((t)/Dt);
			t = t/Dt-j; 
		      // cout<<Dt<<endl;	
			if(j>m+1 || j<0) printf("warning! t=%f out of range!!, Dt = %f, T = %d, m = %d\n",t,Dt,T,m);
			// Calculate square and cube, and pointer to the values to use
			double t2=t*t,t3=t2*t;
			// Calculate the value of the spline function
			if(j==m){bvals[nn] = x[2*j];}
			else{bvals[nn] = x[2*j]*(2*t3-3*t2+1)+x[2*j+1]*(t3-2*t2+t)*Dt+x[2*j+2]*(-2*t3+3*t2)+x[2*j+3]*(t3-t2)*Dt;}
		}
	}
}



#endif

