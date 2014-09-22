function y = baryvalscalar(x,X,Y)

% scalar version--probably so I can just toss this in C++
% at Chebyshev nodes X(j), return Y(j)
N = length(X)-1;
j = find(x==X);
if(j)
    y =Y(j);
else
ds = x-(X);
ds(1) = 2*ds(1);
ds(N+1) = 2*ds(N+1);
y = sum((-1).^[0:N].*Y./(ds))/sum((-1).^[0:N]./(ds));
end