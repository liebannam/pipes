
clear all;
close all;

format short g
%evaluate h(A) at a_j, where
%x_j =-1+c*(a_j)^(2/3)
%x_j are standard Chebyshev nodes 

N = 10;
x = -cos(pi*[0:N]/N); %regular Chebyshev nodes
c = 2*(8/pi)^(2/3);   %coefficient for transform
a = ((1+x)/c).^(3/2); %a_j = A values where we need solution h(A)
xhat = c*a.^(2/3)-1;  %evaluation points 

for(k = 1:N+1)
    th(k) = fzero(@(x) 1/8*(x-sin(x))-a(k), a(k));
end
h = 1/2*(1-cos(th/2));

figure(1)
hold on
plot(a,h,'*')
plot(pi/4-a(N+1:-1:1),1-h(N+1:-1:1),'r*');

%%fit a polynomial on these points
p = polyfit(c*a.^(2/3)-1,h,N);
s = linspace(0,pi/8);
h2 = polyval(p,c*s.^(2/3)-1);
plot(s,h2,'k')
xlabel('A')
ylabel('h')
legend('h(A)', 'odd extension', 'polynomial approx');
title('h(A)')

%%now do it with Chebfun!
f = polyfit(c*a.^(2/3)-1,h,N, domain(-1,1));
fbad = polyfit(a*4/pi-1,h,N, domain(-1,1));
h3 = f(c*s.^(2/3)-1);

%display coefficients
disp('Expansion Coefficients')
disp('c_k denote Chebyshev (from Chebfun), P_k denote basic polynomial (from Matlab polyfit)')
disp('            k  c_k(x)       b_k(x^(2/3))  c_k(x^(2/3)');
fc = chebcoeffs(f);
fcbad = chebcoeffs(fbad);
disp([[0:N]'  fcbad  p([N+1:-1:1])'  fc])



%check out and plot error
for k = 1:length(s)
    u = fzero(@(x) 1/8*(x-sin(x))-s(k), s(k));
    htrue(k) = 0.5*(1-cos(u/2));
    %ds = s(k)^(2/3)-(a.^(2/3));
    %ds(1) = 2*ds(1);
    %ds(N+1) = 2*ds(N+1);
    %errC(k) = htrue(k)- sum((-1).^[0:N].*h./(ds))/sum((-1).^[0:N]./(ds));
    %errD(k) = htrue(k) - baryval(s(k)^(2/3),a.^(2/3), h);  
end
  errA = htrue-h2;
  errB = htrue-h3;
   
figure(2)
plot(s,errA)
hold on
plot(s,errB,'g')
title('errors')
legend('regular old matlab polyfit', 'chebfun')
xlabel('A')
ylabel('fzero-p(A)')

figure(3)
semilogy(N:-1:0, abs(p),'LineWidth', 2)
hold on
semilogy([0:N], abs(fc),'g', 'LineWidth', 2)
title(sprintf('coefficient decay for  N =%d',N))
legend('polyfit', 'chebyshev')
grid on

%timings
disp('Timings comparison for 1e5 evaluations');
s =  linspace(0,pi/8, 1e5)';
a = reshape(a,1,N+1);
h = reshape(h,1,N+1);
disp('My Evaluation')
ss = s.^(2/3);
tic
hs = baryval(s.^(2/3), a.^(2/3), h);
toc
disp('Chebfun Evaluation')
tic
hs2 = f(c*s.^(2/3)-1);
toc
disp('Mean disagreement')
disp(mean(abs(hs2-hs)))

%%for looping--this is really bad now!!
%tic
%hs = zeros(size(s));
% for k =1:length(s)
 %   hs(k) = baryval(s(k).^(2/3), a.^(2/3),h);
 %end
%toc
