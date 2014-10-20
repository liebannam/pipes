
%test scaling to better fit polynomials
close all
clear all
N = 10;
x= -cos(pi*[0:N]/N)';
fy = [];
fys = [];
Y = [];
legs= {};
xr = (x+1)/2;  %regular chebyshev shift
r = linspace(0,1)';
for k = 1:6
    alpha = 2*k/(6);
    y = xr.^(alpha);
    xs = ((x+1)/2).^(1/alpha);
    ys = xs.^(alpha);
    fs = polyfit(x,y,N, domain(-1,1));
    fsb = polyfit(x, ys, N, domain(-1,1));
    fs.coeffs();
    fsb.coeffs();
    Y = [Y y];
    fy = [fy fs(2*r-1)];
    fys = [fys fsb(2*r.^alpha-1)];
    yact = r.^(alpha);
    err = norm(fy(:,k)-yact);
    err2 = norm(fys(:,k)-yact);
    legs{k} = strcat('\alpha = ',sprintf(' %.2f, error = %e error2 = %e', alpha, err,err2));
end


ftest = @(x) sqrt(x);
%ftest = @(x) x.^3-4*x+1;
%ftest = @(x) atan(x*pi);
%test = @(x) tan(x*pi/2);
%ftest = @(x) atan(sin(2*x-1)); %this one is interesting!
%ftest = @(x) sin(pi*x);
%ftest = @(x) cos(4*pi*x);
M = 200;
s = linspace(0,1,200);
alphas = [1:M]'/M
for k =1:length(alphas)
    alpha = alphas(k);
    y = xr.^(alpha);
    xs = ((x+1)/2).^(1/alpha);
    ys = ftest(xs);
    f = polyfit(x, ys, N, domain(-1,1));
    b = polyfit([N:-1:0],log(abs(f.coeffs())'),1);
    err(k) = norm(ftest(s)-f((2*s.^alpha-1)));
    rr(k) = b(1);
end


figure(1)
subplot(1,2,1)
t = linspace(0,1);
plot(t, ftest(t))
xlabel('x')
ylabel('f(x)')
title(func2str(ftest))
subplot(1,2,2)
[ax,p1,p2] = plotyy(alphas,err,alphas,rr,'semilogy','plot');
ylabel(ax(1),'error')
ylabel(ax(2), 'coefficient decay')
grid on
xlabel('\alpha')

figure(2)
plot(xr,Y,'*')
hold on
plot(r,fy)
plot(r,fys, ':')
legend(legs)
title('Fitting x^\alpha')
    