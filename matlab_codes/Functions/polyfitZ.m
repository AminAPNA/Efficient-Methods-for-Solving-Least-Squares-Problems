function [y,H] = polyfitZ(z,f,n,v)
%This function is provided for the Algorithm 1, in the paper
%It Generates an orthonormal basis Q for the Krylov subspace
%It also generates the Hessenberg recurrence matrix H
%It gives the new vector of unknowns y=Q^{H}diag(v)f
%INPUT:
%   z = The vector of order m, including the nodes
%   v = column vector of size m
%   f = the column vector including the data values
Q = v/norm(v);%This constructs Q and H using Arnoldi iteration
H = zeros(n+1,n);
for k = 1:n
    q = z.*Q(:,k);
    for j = 1:k
        H(j,k) = Q(:,j)'*q;
        q = q - H(j,k)*Q(:,j);
    end
    H(k+1,k) = norm(q);
    for j = 1:k  %Reorthogonalization loop
        mu=Q(:,j)'*q;
        q = q - mu*Q(:,j);
        H(j,k) = H(j,k)+mu;
    end
    H(k+1,k) = norm(q);
    Q=[Q q/H(k+1,k)];
end
y = Q\(v.*f);%It computes y using the matrix Q obtained by Arnoldi iteration
end
