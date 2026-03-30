function [y,H,K] = ratfitZ(z,Pole,f,n,v)
%This function is provided for the Algorithm 5, in the paper
%It Generates an orthonormal basis Q for the rational Krylov subspace
%It also generates the Hessenberg pencil (H,K)
%It gives the new vector of unknowns y=Q^{H}diag(v)f
%INPUT:
%   z = The vector of order m, including the nodes
%   Pole = The column vector including n poles
%   v = column vector of size m
%   f = the column vector including the data values
m = length(z);%This constructs Q and H and K using rational Arnoldi iteration
Q = v/norm(v);
H = zeros(n+1,n);
for k = 1:n
    if k == 1
        q = (ones(m,1)./(z-Pole(1))).*Q(:,1);
    else
        q = ((z-Pole(k-1))./(z-Pole(k))).*Q(:,k);
    end
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
E = eye(n+1,n); E(1,1) = 0;%As the first shift is (-\infty), we use the fraction form of the first shift as (-1/0)
K = H - E;
H = H*diag(Pole) - eye(n+1,n)*diag([-1,Pole(1:n-1)]);
y = Q\(v.*f);%It computes y using the matrix Q obtained by rational Arnoldi iteration
end