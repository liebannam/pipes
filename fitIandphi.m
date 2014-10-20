clear all
close all

tic 
N = 10;
Af = pi/4;
eps = Af*1e-5;
g = 9.8;
M = 3*2^2;

tx = @(x,p,ym) ym*((x+1)/2).^(1/p);
ty = @(y,p,ym) 2*(y/ym).^(p)-1;

alphas = [(2:M)/M];
%alphas = [1/3];
x = -cos(pi*[0:N]/N)'; %regular Chebyshev nodes
ca = 2*(8/pi)^(2/3);   %coefficient for transform for A
a = ((1+x)/ca).^(3/2); %a_j = A values where we need solution for h(A) interpolation 
ch = 2^(3/2);            %coefficient for transform for h

for(k = 1:N+1)
    th(k) = fzero(@(x) 1/8*(x-sin(x))-a(k), a(k));
end
h = 1/2*(1-cos(th/2))';
H = ((1+x)/ch).^(2);  
T = 2*acos(1-2*H);
A = 1/8*(T-sin(T));

f = polyfit(2*(a/(pi/8)).^(2/3)-1,h,N, domain(-1,1));
finv = polyfit(2*(H/.5).^(1/2)-1, A, N, domain(-1,1));

D = 1;
hofA = @(x) f(ca*x.^(2/3)-1).*(x<=pi/8) + (1-f(ca*(pi/4-x).^(2/3)-1)).*(x>pi/8);
Aofh = @(x) finv(ch*x.^(1/2)-1).*(x<=0.5)+(pi/4-finv(ch*(1-x).^(1/2)-1)).*(x>0.5);
l = @(x) sqrt(1.-(1-2*x).^2);
fphi = @(x) sqrt(1./(x.*l(hofA(x))));
fI = @(x,A) (hofA(A)-x).*(l(x));
Ie = @(y) g/12.*((3*D*D-4*D*y+4*y.*y).*sqrt(y.*(D-y))-3*D*D*(D-2*y).*atan(sqrt(y)./sqrt(D-y)));
powphi = @(t) sqrt(1*3.*D/8.)*(t-1/80.*t.^3+19./448000.*t.^5+1./10035200.*t.^7+491./(27.*70647808000.)*t.^9);

for j = 1:length(alphas)
    p = alphas(j);
    a1 = pi/8*((x+1)/2).^(1/p);
    a2 = -pi/8*((x+1)/2).^(1/p)+pi/4;

    ta1 = @(a) 2*(8*a/pi).^(p)-1;
    ta2 = @(a) 2*((pi/4-a)*8/pi).^(p)-1;

    phimax = quad(fphi, 0,pi/4,1e-12);
    for k = 1:N+1 
        I1(k) = Ie(hofA(a1(k)));
        I2(k) = Ie(hofA(a2(k)));
        phi1(k) = quad(fphi, 0,a1(k),1e-12);
        phi2(k) = quad(fphi,0,a2(k),1e-12);
    end
    I1(1) = 0;
    oops = [find(isnan(phi1)) find(isinf(phi1))];
    phi1(oops) = fphi(a(oops)).*a(oops);
    phi1(1) =0;
    fI1 = polyfit(ta1(a1), I1', N, domain(-1,1));
    fI2 = polyfit(ta2(a2), I2', N,domain(-1,1));    
    fp1 = polyfit(ta1(a1), phi1', N, domain(-1,1));
    fp2 = polyfit(ta2(a2), phi2',N,domain(-1,1));
    IofA = @(x) fI1(ta1(x)).*(x<=pi/8) + (fI2(ta2(x))).*(x>pi/8);
    phiofA = @(x) fp1(ta1(x)).*(x<=pi/8)+(fp2(ta2(x))).*(x>pi/8);
    phimax1 = phiofA(pi/8);
    dphi = phimax-phimax1;
    
    tx5 = @(x) phimax1*((x+1)/2);
    ta5 = @(a) 2*(a/phimax1)-1;
%    tx5 = @(x) phimax1*((x+1)/2).^(1/p);
%    ta5 = @(a) 2*(a/phimax1).^(p)-1;
    tx6 = @(x) -dphi*((x+1)/2).^(1/p)+phimax;
    ta6 = @(a) 2*((phimax-a)/dphi).^(p)-1;
    
    phis1 = tx5(x);
    phis2 = tx6(x);
    for k = 2:N
        as1(k) = fzero(@(x) phiofA(x)-phis1(k),[0,pi/8]);
        as2(k) = fzero(@(x) phiofA(x)-phis2(k),[pi/8,pi/4]);
    end
    as1(1) = 0;
    as1(N+1) = pi/8;
    as2(1+N) = pi/8;
    as2(1) = pi/4;

    fphiinv1 = polyfit(ta5(phis1),(as1.^(p))', N,domain(-1,1));
    fphiinv1 = polyfit(ta5(phis1),(as1)', N,domain(-1,1));
    fphiinv2 = polyfit(ta6(phis2),(as2)', N,domain(-1,1));
    Aofphi = @(x) (fphiinv1(ta5(x))).*(x<phimax1)+fphiinv2(ta6(x)).*(x>=phimax1);

    t1 = linspace(0,pi,250);
    t2 = linspace(pi,2*pi,250);
    pt1 = linspace(0,phimax1,250);
    pt2 = linspace(phimax1,phimax,250);
    
    at1 = 1/8*(t1-sin(t1));
    at2 = 1/8*(t2-sin(t2));
    %ht = 1/2*(1-cos(t/2));
    
    errs(j,1) = norm(Ie(hofA(at1))-IofA(at1));
    errs(j,2) = norm(Ie(hofA(at2))-IofA(at2));
    errs(j,3) = norm(phiofA(Aofphi(pt1))-pt1);
    errs(j,4) = norm(Aofphi(phiofA(at1))-at1);
    errs(j,5) = norm(phiofA(Aofphi(pt2))-pt2);
    errs(j,6) = norm(Aofphi(phiofA(at2))-at2);
    %b1 = polyfit([N:-1:0],log(abs(f.coeffs())'),1);
    %b2 = polyfit([N:-1:0],log(abs(finv.coeffs())'),1);
   % r(j,1) = b1(1);
   % r(j,2) = b1(2);
   fci1 = fI1.coeffs();
   fci2 = fI2.coeffs();
   fcp1 = fp1.coeffs();
   fcp2 = fp2.coeffs();
   fcpi1 = fphiinv1.coeffs();
   fcpi2 = fphiinv2.coeffs();
   r(j,1) = abs(fci1(1))
   r(j,2) = abs(fci2(1))
   r(j,3) = abs(fcp1(1))
   r(j,4) = abs(fcp2(1))
   r(j,5) = abs(fcpi1(1));
   r(j,6) = abs(fcpi2(1));
end
figure(1)
subplot(2,1,1)
h=semilogy(alphas,errs)
set(h, 'LineWidth', 2)
set(gca, 'YTick', [1e-12,1e-10,1e-8,1e-6,1e-4,1e-2])
set(gca, 'XTick', [1/6,1/4,1/3,1/2,2/3,3/4,5/6,1])
xlab = {'1/6','1/4', '1/3', '1/2', '2/3','3/4','5/6', '1  '};
set(gca, 'XTickLabel',xlab)
grid on
title('Fitting with scaling x^\alpha')
legend('I_1(A)','I_2(A)','\phi_1(A(\phi_1))-\phi_1','A(\phi_1(A))-A ','\phi_2(A(\phi_2))-\phi_2', 'A(\phi_2(A))-A')
xlabel('\alpha')
ylabel('Error')

subplot(2,1,2)
%plot(alphas, r)
h = semilogy(alphas,r)
set(h, 'LineWidth', 2)
set(gca, 'XTick', [1/6,1/4,1/3,1/2,2/3,3/4,5/6,1])
xlab = {'1/6','1/4', '1/3', '1/2', '2/3','3/4','5/6', '1  '};
set(gca, 'XTickLabel',xlab)
grid on
legend('I_1(A)','I_2(A)','\phi_1(A)','\phi_2(A)','\phi^{-1}_1', '\phi^{-1}_2')
xlabel('\alpha')
ylabel('Coefficient Decay')
toc