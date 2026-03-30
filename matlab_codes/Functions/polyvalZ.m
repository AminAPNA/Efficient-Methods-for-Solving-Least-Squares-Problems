function p = polyvalZ(y,H,x,p0)
%This function is provided for the Algorithm 2, in the paper
%It Generates the least squares solution at an arbitrary set of sampling nodes (column vector x of order M)
%INPUT:
%   y = the column vector y computed through the function polyfitZ (Algorithm 1)
%   H = Hesssenberg matrix obtained from the function polyfitZ
%   x = the column vector x of order M including the sampling nodes
%   p0 = the orthonormal polynomial of degree 0 directly computed through the inner product
M = length(x);
U = zeros(M,1);
n = size(H,2);
U(:,1) = p0;%U is the Mx(n+1) matrix including the orthonormal polynomials at sampling nodes
for k = 1:n
    u = x.*U(:,k);
    for j = 1:k
        u = u - H(j,k)*U(:,j);
    end
    U = [U u/H(k+1,k)];
end
p = U*y;%The vector including p(x_{j})
end
