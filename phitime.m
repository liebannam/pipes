
clear all
close all

N = 15;
Af = pi/4;
eps = Af*1e-5;
g = 9.8;

x = -cos(pi*[0:N]/N)'; %regular Chebyshev nodes
ca = 2*(8/pi)^(2/3);   %coefficient for transform for A
a = ((1+x)/ca).^(3/2); %a_j = A values where we need solution for h(A) interpolation 
ch = 2^(3/2);            %coefficient for transform for h
H = ((1+x)/ch).^(2);   %h_j = h values where we need solution ofr A(h) interpolation


for(k = 1:N+1)
    th(k) = fzero(@(x) 1/8*(x-sin(x))-a(k), a(k));
end
h = 1/2*(1-cos(th/2))';
T = 2*acos(1-2*H);
A = 1/8*(T-sin(T));


%%now do it with Chebfun!
f = polyfit(ca*a.^(2/3)-1,h,N, domain(-1,1));
finv = polyfit(ch*H.^(1/2)-1, A, N, domain(-1,1))

D = 1;
hofA = @(x) f(ca*x.^(2/3)-1).*(x<=pi/8) + (1-f(ca*(pi/4-x).^(2/3)-1)).*(x>pi/8);
Aofh = @(x) finv(ch*x.^(1/2)-1).*(x<=0.5)+(pi/4-finv(ch*(1-x).^(1/2)-1)).*(x>0.5);
l = @(x) sqrt(1.-(1-2*x).^2);
fphi = @(x) sqrt(1./(x.*l(hofA(x))));
fI = @(x,A) (hofA(A)-x).*(l(x));
Ie = @(y) g/12.*((3*D*D-4*D*y+4*y.*y).*sqrt(y.*(D-y))-3*D*D*(D-2*y).*atan(sqrt(y)./sqrt(D-y)));

alpha1 = 1/3;
alpha2 = 1/3;
alpha3 = 1/3;
alpha4 = 1/3;
%a =  pi/4*((x+1)/2).^(1/alpha);
tx1 = @(x) pi/8*((x+1)/2).^(1/alpha1);
ta1 = @(a) 2*(8*a/pi).^(alpha1)-1;
tx2 = @(x) -pi/8*((x+1)/2).^(1/alpha2)+pi/4;
ta2 = @(a) 2*((pi/4-a)*8/pi).^(alpha2)-1


phimax = quad(fphi, 0,pi/4);

tx3 = @(x) pi/8*((x+1)/2).^(1/alpha3);
ta3 = @(a) 2*(8*a/pi).^(alpha3)-1;
%tx4 = @(x) pi/8*((x+1)/2).^(1/alpha4)+pi/8;
%ta4 = @(a) 2*((a-pi/8)*8/pi).^(alpha4)-1;

tx4 = @(x) -pi/8*((x+1)/2).^(1/alpha4)+pi/4;
ta4 = @(a) 2*((pi/4-a)*8/pi).^(alpha4)-1


%ta = @(a) 16/pi*(a-pi/8)-1;
%tx = @(x) pi*(x+1)/16+pi/8;

%tx = @(x) asinh((x+1)/2*sinh(2))*pi/8;
%ta = @(a) 2*sinh(a*8/pi)/sinh(2)-1;
%tx = @(x) log(((exp(pi/4)-1)*x+1+exp(pi/4))/2);
%ta = @(a) (2*exp(a)-exp(pi/4)-1)/(exp(pi/4)-1);
 %tx = @(x) 1/2*(cosh((x+1)/2*acosh(pi/2+1))-1);
 %ta = @(a) 2*(acosh(2*a+1))/acosh(pi/2+1)+1
a = tx1(x);
a2 = tx2(x);
a3 = tx3(x);
a4 = tx4(x);


phi(1) = 0;
I(1) = 0;
h2 = linspace(0,1);
for k = 1:length(a)
    phi(k) = quad(fphi, 0,a3(k));
    phi2(k) = quad(fphi,0,a4(k));
    y = hofA(a(k));
    I(k) = Ie(y);
    I2(k) = Ie(hofA(a2(k)));
    %I2(k) = g*quad(@(t)fI(t,a(k)), 0, hofA(a(k)),1e-12);
    %a2(k) = quad(@(x)sqrt(1-(1-2*x).^2), 0, h2(k),1e-16);
end
oops = find(isnan(phi)) ==1;
phi(oops) = fphi(oops).*a(oops);
%fI = polyfit(4/pi*2*a-1, I', N, domain(-1,1))
fI1 = polyfit(ta1(a), I', N, domain(-1,1));
fI2 = polyfit(ta2(a2), I2', N,domain(-1,1));
fp = polyfit(ta3(a3), phi', N, domain(-1,1));
fp2 = polyfit(ta4(a4), phi2',N,domain(-1,1));
q = rand(1,250)*pi/8;
q2 = rand(1,250)*pi/8+pi/8;
p = rand(1,250);

t = 2*acos(1-2*p);
plot(q,Ie(hofA(q))-fI1(ta1(q)),'*')
hold on
plot(q2, Ie(hofA(q2))-fI2(ta2(q2)),'g*')
disp('       f              finv         fI1       fI2       fp     fp2 ')
disp([f.coeffs()  finv.coeffs() fI1.coeffs() fI2.coeffs() fp.coeffs()  fp2.coeffs()])
disp('errors')
e1 = norm(Aofh(p)-1/8*(t-sin(t)));
e2 = norm(hofA(1/8*(t-sin(t)))-p);
e3 = norm(Ie(hofA(q))-fI1(ta1(q)));
e4 = norm(Ie(hofA(q2)) - fI2(ta2(q2)));
powphi = @(t) sqrt(1*3.*D/8.)*(t-1/80.*t.^3+19./448000.*t.^5+1./10035200.*t.^7+491./(27.*70647808000.)*t.^9);
phiofA = @(x) fp(ta3(x)).*(x<pi/8) +fp2(ta4(x)).*(x>=pi/8);
IofA = @(x) fI1(ta1(x)).*(x<pi/8)+fI2(ta2(x)).*(x>=pi/8);

disp([e1   e2   e3   e4])
s = linspace(0,pi/8,250);
%s = linspace(pi/8,pi/4,250);
% figure(2)
% plot(s,fI1(ta1(s)))
% hold on
% plot(s+pi/8, fI2(ta2(s+pi/8)),'g')
% title('I(A)')
% figure(3)
% plot(s,fp(ta3(s)))
% hold on
% plot(s+pi/8, fp2(ta4(s+pi/8)),'g')
% t = linspace(0,2*pi);
% plot(1/8*(t-sin(t)),powphi(t),'r')
% title('\phi(A)')

alpha5 = 1/2;
alpha6 = 1/3;
phimax1 = phiofA(pi/8);
dphi = phimax-phimax1;
tx5 = @(x) phimax1*((x+1)/2).^(1/alpha5);
ta5 = @(a) 2*(a/phimax1).^(alpha5)-1;

tx6 = @(x) -dphi*((x+1)/2).^(1/alpha6)+phimax;
ta6 = @(a) 2*((phimax-a)/dphi).^(alpha6)-1

phis1 = tx5(x);
phis2 = tx6(x);
for k = 2:N
    as1(k) = fzero(@(x) phiofA(x)-phis1(k),[0,pi/4]);
    as2(k) = fzero(@(x) phiofA(x)-phis2(k),[0,pi/4]);
end
as1(1) = 0;
as1(N+1) = pi/8;
as2(N+1) = pi/8;
as2(1) = pi/4;

fphiinv1 = polyfit(ta5(phis1),as1', N,domain(-1,1));
fphiinv2 = polyfit(ta6(phis2),as2', N,domain(-1,1));
[fphiinv1.coeffs()  fphiinv2.coeffs()]


Aofphi = @(x) fphiinv1(ta5(x)).*(x<phimax1)+fphiinv2(ta6(x)).*(x>=phimax1);
figure(2) 
A = linspace(0,pi/4,500);
plot(A/pi*4, hofA(A), 'g', 'LineWidth', 2)
hold on
plot(A/pi*4, IofA(A), 'b','LineWidth', 2)
plot(A/pi*4,phiofA(A),'k','LineWidth', 2)
legend('h', 'I', 'phi');
xlabel('A/Af')
set(gca, 'XTick', [0:.25:1])
set(gca, 'YTick', [0:.5:4])
grid on
axis([-.05, 1.05,0,4])

disp('scalings')
disp('h(A) ~sum_k ak*Tk(A^(2/3))')
disp('A(h) ~sum_k bk*Tk(h^(1/2))')
disp(sprintf('A<pi/8: I1(A) ~sum_k bk*Tk(A^(1/%d))',1/alpha1))
disp(sprintf('A>= pi/8: I2(A) ~sum_k bk*Tk(A^(1/%d))',1/alpha2))
disp(sprintf('A<pi/8: phi1(A) ~sum_k bk*Tk(A^(1/%d))',1/alpha3))
disp(sprintf('A>= pi/8: phi2(A) ~sum_k bk*Tk(A^(1/%d))',1/alpha4))
t = linspace(0,phimax,1000)';
aa = linspace(0,pi/4,1000)';
disp('errors in phi(A(phi))etc')
[(phiofA(Aofphi(t))-t)  (Aofphi(phiofA((aa)))-aa)]
disp('errors in hofA(Aofh(h)) etc')
hh = linspace(0,1)';
close all
aha = Aofh(hofA(aa))-aa;
hah = (hofA(Aofh(hh))-hh);
pap = (phiofA(Aofphi(t)))-t;
apa = Aofphi(phiofA(aa))-aa;
oops(1) = norm(imag(aha));
oops(2) = norm(imag(aha));
oops(3) = norm(imag(pap));
oops(4) = norm(imag(apa));
semilogy(aa/max(aa),abs(aha));
hold on
semilogy(hh/max(hh),abs(hah),'g');
semilogy(t/max(t),abs(pap),'k');
semilogy(aa/max(aa), abs(apa),'m');
legs = {'h(A(h))-h', 'A(h(A))-A','\phi(A(\phi))-\phi', 'A(\phi(A))-A'};
legend(legs)
disp('imaginary bits?')
disp(legs)
disp(oops)
