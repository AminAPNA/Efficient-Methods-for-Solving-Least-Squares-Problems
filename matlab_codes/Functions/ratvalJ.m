function r = ratvalJ(y,H,K,X,S,r0)
%This function is provided for the Algorithm 8, in the paper
%It Generates the least squares solution at an arbitrary set of sampling nodes (column vector x of order \tau)
%INPUT:
%   y = the column vector y computed through the function ratfitJ (Algorithm 7)
%   H = Hesssenberg matrix obtained from the function ratfitJ
%   K = Hesssenberg matrix obtained from the function ratfitJ
%   X = the Jordan-like matrix constructed by the sampling nodes x
%   S = The vector of order \tau including the sampling derivatives \tilde{s}_{j}
%   r0 = the orthonormal rational function of degree 0 directly computed through the inner product
tau = length(S);
M = size(X,1);
n = size(H,2);
U = [];%U is the Mx(n+1) matrix including the derivatives up to order \tilde{s}_{j} of the Sobolev orthonormal rational functions at sampling nodes
for j = 1:tau
    U1 = zeros(S(j)+1,1);
    U1(S(j)+1,1) = r0;
    U = [U;U1];
end
for k = 1:n
    u = 0;
    for j = 1:k
        u = u + K(j,k)*X*U(:,j) - H(j,k)*U(:,j);
    end
    U = [U (H(k+1,k)*eye(M)-K(k+1,k)*X)\u];
end
r = U*y;%The vector including the approximations
end