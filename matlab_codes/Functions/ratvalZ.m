function r = ratvalZ(y,H,K,x,r0)
%This function is provided for the Algorithm 6, in the paper
%It Generates the least squares solution at an arbitrary set of sampling nodes (column vector x of order M)
%INPUT:
%   y = the column vector y computed through the function ratfitZ (Algorithm 5)
%   H = Hesssenberg matrix obtained from the function ratfitZ
%   K = Hesssenberg matrix obtained from the function ratfitZ
%   x = the column vector x of order M including the sampling nodes
%   r0 = the orthonormal rational function of degree 0 directly computed through the inner product
M = length(x);
n = size(H,2);
U = zeros(M,1);
U(:,1) = r0;%U is the Mx(n+1) matrix including the orthonormal rational functions at sampling nodes
for k = 1:n
    u = 0;
    for j = 1:k
        u = u + H(j,k)*U(:,j) - K(j,k)*x.*U(:,j);
    end
    U = [U u./(K(k+1,k)*x-H(k+1,k))];
end
r = U*y;%The vector including r(x_{j})
end