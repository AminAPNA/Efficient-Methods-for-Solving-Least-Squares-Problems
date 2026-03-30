function [y,H] = polyfitJ(J,f,n,w,v)
%This function is provided for the Algorithm 3, in the paper
%It Generates an orthonormal basis Q for the Krylov subspace
%It also generates the Hessenberg recurrence matrix H
%It gives the new vector of unknowns y=Q^{H}diag(w)f
%INPUT:
%   J = mxm matrix
%   w = the column vector of order m including weights and factors compensating for the repeated differentiation which introduces factorial terms
%   v = column vector of size m
%   f = the column vector including the data points
Q = v/norm(v);%This constructs Q and H using Arnoldi iteration
H = zeros(n+1,n);
for k = 1:n
    q = J*Q(:,k);
    for j = 1:k
        H(j,k) = Q(:,j)'*(q);
        q = q - H(j,k)*Q(:,j);
    end
    H(k+1,k) = norm(q);
    for j = 1:k  %Reorthogonalization loop
        mu=Q(:,j)'*q;
        q = q - mu*Q(:,j);
        H(j,k) = H(j,k)+mu;
    end
    H(k+1,k) = norm(q);
    Q = [Q q/H(k+1,k)];
end
y = Q\(w.*f);%It computes y using the matrix Q obtained by Arnoldi iteration
end