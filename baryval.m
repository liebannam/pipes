function y = baryval(x,X,Y)

% Barycentric Evaluation for Chebyshev interpolant of y values Y and at Chebyshev nodes X 
% vectorized verson
% evaluate at points in vector x
% scale invariant (eg any chebyshev nodes on any interval will work)
% but it will be terrible if x is outside of bounds of X

% following Trefethen's Approximation Theory, pg 35
% note: about a factor of three slower than the chebfun implementation :(

% x is reshaped to be an Mx1 column vector
% X,Y are reshaped to be 1xN+1 row vectors

%kind of whack in testing for Chebyshv points but seems to work
N = length(X)-1;
M = length(x);
where = zeros(M,1);
for(i=1:N+1)
    f = find(x==X(i));
    if f
        where(f) =i;
    end
end

x = reshape(x,M,1);
X = reshape(X,1,N+1);
Y = reshape(Y,1,N+1);
ds = kron(x,ones(1,N+1))-(kron(X,ones(M,1)));
v = (-1).^(0:N);
v(1) = 2*v(1);
v(N+1)= 2*v(N+1); 
ds = ds*diag(v);
y = sum(kron(Y, ones(M,1))./(ds),2)./sum(1./(ds),2);

%%if there were any Chebyshev points in x, y has NaNs, need to replace with
%%correct values

if(sum(where)>0)    
    y(where>0) = Y(where(where>0))';
end


