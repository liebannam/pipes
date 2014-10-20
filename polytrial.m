%%%eeeee it works so much better when you scale it right!!!!
%%correct scalings for this problem:

%evaluate h(A) at a_j, where
%x_j =-1+ca*(a_j)^(2/3)

%evaluate A(h) at h_j, where
%x_j = -1+ch*h_j^(1/2)

%where x_j are standard Chebyshev nodes 

%%odd extentions:
%hofA should be the following function:
%hofA(A) = fa(A)  (A<=pi/8)
%         =1-fa(pi/4-A)   (pi/8<A<=pi/4)

%likewise Aofh is
%Aofh(h) = fh(h)    (h<=1/2)
%        = pi/4-fh(1-h)  (1/2<h<=1)

%%Still need: c, phi, eta (sigh/facepalm)
%c scaling seems to by 1/3 but the singular right endpoint is seriously
%cramping my style 
%eta scaling is 1/3 for a in (0,pi/8) but it looks bad for (0,pi/4)...
%transition height strategy: find At such that Aofh(ht) = At (is in theory
%farther from Af than you get from using 1/8(theta-sin(theta))
clear all;
close all;

format short g

N = 15;
Af = pi/4;
eps = Af*1e-5;

x = -cos(pi*[0:N]/N)'; %regular Chebyshev nodes
ca = 2*(8/pi)^(2/3);   %coefficient for transform for A
a = ((1+x)/ca).^(3/2); %a_j = A values where we need solution for h(A) interpolation 
ch = 2^(3/2);            %coefficient for transform for h
H = ((1+x)/ch).^(2);   %h_j = h values where we need solution ofr A(h) interpolation
cc = 2*(8/pi).^(1/3);
ac = ((1+x)/cc).^(3);
%ac2 = ((1+x)/cc2).^(3); %c evaluation points for A>pi/8

a2 = (1+x)*pi/16;  %regular linear transform--will be terrible!
H2 = (1+x)/4;

%c2 = 2^(5/2);
%H = ((1+x)/c2).^(2/3);  %h values to use for H(A) interpolation(?)
%H = ((1+x)/k).^2;

%xhat = c*a.^(2/3)-1;  %evaluation points 

for(k = 1:N+1)
    th(k) = fzero(@(x) 1/8*(x-sin(x))-a(k), a(k));
    th2(k) = fzero(@(x) 1/8*(x-sin(x))-a2(k), a2(k));
    th3(k) = fzero(@(x) 1/8*(x-sin(x))-ac(k), ac(k));
    %th4(k) = fzero(@(x) 1/8*(x-sin(x))-ac2(k), ac2(k));
end
h = 1/2*(1-cos(th/2))';
T = 2*acos(1-2*H);
T2 = 2*acos(1-2*H2);
A = 1/8*(T-sin(T));
A2 = 1/8*(T2-sin(T2));
h22 = 1/2*(1-cos(th2/2))';
c = sqrt(9.8*ac./sin(th3/2)');
c(1) = 0;

figure(1)
hold on
plot(a,h,'*')
plot(pi/4-a(N+1:-1:1),1-h(N+1:-1:1),'r*');

%%fit a polynomial on these points
p = polyfit(ca*a.^(2/3)-1,h,N);
s = linspace(0,pi/8);
h2 = polyval(p,ca*s.^(2/3)-1);
plot(s,h2,'k')

xlabel('A')
ylabel('h')
legend('h(A)', 'odd extension', 'polynomial approx');
title('h(A)')

%%now do it with Chebfun!
f = polyfit(ca*a.^(2/3)-1,h,N, domain(-1,1));
fbad = polyfit(a2,h22,N, domain(0,pi/8));
%finv = polyfit(c2*H.^(3/2)-1, A, N, domain(-1,1));
%finv = polyfit(H, A, N, domain(0,.5));
finv = polyfit(ch*H.^(1/2)-1, A, N, domain(-1,1));
finvbad  = polyfit(H2,A2,N, domain(0,.5)); 
h3 = f(ca*s.^(2/3)-1);

%display coefficients
disp('Expansion Coefficients')
disp('c_k denote Chebyshev (from Chebfun), P_k denote basic polynomial (from Matlab polyfit)')
disp('            k  c_k(x)       b_k(x^(2/3))  c_k(x^(2/3)');
%fc = chebcoeffs(f);
%fcbad = chebcoeffs(fbad);
fc = f.coeffs();
fcbad = fbad.coeffs();
fic = finv.coeffs();
fc = fc(N+1:-1:1);
fcbad = fcbad(N+1:-1:1);
fic = fic(N+1:-1:1);
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
  errA = abs(htrue-h2);
  errB = abs(htrue-h3);
  errC = abs(htrue-fbad(s));
  
figure(2)
semilogy(s,errA)
hold on
semilogy(s,errB,'g')
title('errors')
semilogy(s,errC,'r')
legend('regular old matlab polyfit', 'chebfun', 'not scaled')
xlabel('A')
ylabel('fzero-p(A)')

figure(3)
semilogy(N:-1:0, abs(p),'LineWidth', 2)
hold on
semilogy([0:N], abs(fc),'g', 'LineWidth', 2)
title(sprintf('coefficient decay for  N =%d',N))
semilogy([0:N], abs(fic), 'm')
legend('h(A) polyfit', 'h(A) chebyshev', 'A(h) chebyshev')
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
hs2 = f(ca*s.^(2/3)-1);
toc
disp('Mean disagreement')
disp(mean(abs(hs2-hs)))
thetas = 2*(acos(1-2*hs));
disp('Mean Error from back-solving for A')
disp(norm(abs(1/8*(thetas-sin(thetas))-s)))

%hofA = @(x) f(ca*x.^(2/3)-1);
%Aofh = @(x) finv(ch*x.^(1/2)-1);
hofA = @(x) f(ca*x.^(2/3)-1).*(x<=pi/8) + (1-f(ca*(pi/4-x).^(2/3)-1)).*(x>pi/8);
Aofh = @(x) finv(ch*x.^(1/2)-1).*(x<=0.5)+(pi/4-finv(ch*(1-x).^(1/2)-1)).*(x>0.5);
figure(1)
t2 = linspace(0,pi/4);
plot(t2,hofA(t2),'m');

t =linspace(0,pi);
atrue = 1/8*(t-sin(t));
htrue = 1/2*(1-cos(t/2));
norm(Aofh(htrue)-atrue);
norm(hofA(atrue)-htrue);
    

% r = linspace(0,1)';
% fy = [];
% Y = [];
% legs= {};
% xr = (x+1)/2;  %regular chebyshev shift
% for k = 1:6
%     alpha = 2*k/(6);
%     y = xr.^(alpha);
%     fs = polyfit(xr,y,N, domain(-1,1));
%     fsb = polyfit(x, y, N, domain(-1,1));
%     fs.coeffs();
% `    %fsb.coeffs();
%     Y = [Y y];
%     fy = [fy fs(r)];
%     err = norm(fs(r)-abs(r).^(alpha).*sign(r));
%     legs{k} = strcat('\alpha = ',sprintf(' %.2f, error = %.4f', alpha, err));
% end
% figure(4)
% plot(x,Y,'*')
% hold on
% plot(r,fy)
% legend(legs)
%title('Fitting x^\alpha')    