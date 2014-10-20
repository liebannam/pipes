clear all
close all

N = 15;
Af = pi/4;
eps = Af*1e-5;
g = 9.8;
M = 3*2^6;

tx = @(x,p,ym) ym*((x+1)/2).^(1/p);
ty = @(y,p,ym) 2*(y/ym).^(p)-1;
alphas = (1:M)/M;
x = -cos(pi*[0:N]/N)'; %regular Chebyshev nodes
for j = 1:length(alphas)
    p = alphas(j);
    a = tx(x,p,pi/8); %a_j = A values where we need solution for h(A) interpolation 
    H = tx(x,p,1/2);   %h_j = h values where we need solution ofr A(h) interpolation
    ta = @(a) ty(a,p,pi/8);
    th = @(h) ty(h,p,1/2);
    for k = 1:N+1
        the(k) = fzero(@(x) 1/8*(x-sin(x))-a(k), a(k));
    end
    h = 1/2*(1-cos(the/2))';
    T = 2*acos(1-2*H);
    A = 1/8*(T-sin(T));
    
    f = polyfit(ta(a),h,N, domain(-1,1));
    finv = polyfit(th(H), A, N, domain(-1,1));
   % p
   % f.coeffs()
   % finv.coeffs()
    hofA = @(x) f(ta(x)).*(x<=pi/8) + (1-f(ta(pi/4-x))).*(x>pi/8);
    Aofh = @(x) finv(th(x)).*(x<=0.5)+(pi/4-finv(th(1-x))).*(x>0.5);
    t = rand(500,1)*2*pi;
    t = linspace(0,2*pi,500);
    at = 1/8*(t-sin(t));
    ht = 1/2*(1-cos(t/2));
    errs(j,1) = norm(ht-hofA(at));
    errs(j,2) = norm(at-Aofh(ht));
    b1 = polyfit([N:-1:0],log(abs(f.coeffs())'),1);
    b2 = polyfit([N:-1:0],log(abs(finv.coeffs())'),1);
   % r(j,1) = b1(1);
   % r(j,2) = b1(2);
   fc = f.coeffs();
   fci = finv.coeffs();
   r(j,1) = abs(fc(1));
   r(j,2) = abs(fci(1));
   
end
figure(1)

subplot(2,1,1)
h=semilogy(alphas,errs)
set(h, 'LineWidth', 2)
set(gca, 'XTick', [1/6,1/4,1/3,1/2,2/3,1])
set(gca, 'YTick', [1e-12,1e-10,1e-8,1e-6,1e-4,1e-2])
xlab = {'1/6','1/4', '1/3', '1/2', '2/3', '1  '}
set(gca, 'XTickLabel',xlab)
grid on
title('Fitting with scaling x^\alpha')
legend( 'h(A)', 'A(h)')
xlabel('\alpha')
ylabel('Error')

subplot(2,1,2)
%plot(alphas, r)
h = semilogy(alphas,r)
set(h, 'LineWidth', 2)
set(gca, 'XTick', [1/6,1/4,1/3,1/2,2/3,1])
xlab = {'1/6','1/4', '1/3', '1/2', '2/3', '1  '}
set(gca, 'XTickLabel',xlab)
grid on
legend('h(A)', 'A(h)')
xlabel('\alpha')
ylabel('Coefficient Decay')
