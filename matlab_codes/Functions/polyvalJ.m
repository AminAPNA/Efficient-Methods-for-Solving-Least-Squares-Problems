function p = polyvalJ(y,H,X,S,p0)
%This function is provided for the Algorithm 4, in the paper
%It Generates the least squares solution at an arbitrary set of sampling nodes (column vector x of order \tau)
%INPUT:
%   y = the column vector y computed through the function polyfitJ (Algorithm 3)
%   H = Hesssenberg matrix obtained from the function polyfitJ
%   X = the Jordan-like matrix constructed by the sampling nodes x
%   S = The vector of order \tau including the sampling derivatives \tilde{s}_{j}
%   p0 = the orthonormal polynomial of degree 0 directly computed through the inner product
tau = length(S);
n = size(H,2);
U = [];%U is the Mx(n+1) matrix including the derivatives up to order \tilde{s}_{j} of the Sobolev orthonormal polynomials at sampling nodes
for j = 1:tau
    U1 = zeros(S(j)+1,1);
    U1(S(j)+1,1) = p0;
    U = [U;U1];
end
for k = 1:n
    u = X*U(:,k);
    for j = 1:k
        u = u - H(j,k)*U(:,j);
    end
    U = [U u/H(k+1,k)];
end
p = U*y;%The vector including the approximations
end