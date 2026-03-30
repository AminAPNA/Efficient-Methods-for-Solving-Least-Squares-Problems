%This code is done for Example 4, in which a Sobolev rational least squares problem is considered
%It Generates two plots related to Example 4 (Figure 4)
%INPUT:
%   J = M1xM1 matrix (We consider M1 nodes)
%   v = the column vector of order M1 including weights and factors compensating for the repeated differentiation which introduces factorial terms
%   v2 = column vector of size M1
%   ff1 = the column vector including the data points
%    t = the vector of poles
%   S = the Jordan-like matrix constructed by the sampling nodes s1
%   s2 = The vector of order m including the sampling derivatives \tilde{s}_{j}
%   r0 = the Sobolev orthonormal rational function of degree 0 directly computed through the inner product


%First we compute the matrix Q, H, K and also vector y through the function ratfitJ
%Then by taking a set of sampling nodes, we evaluate the approximations by ratvalJ.
%At last we compute the errors at sampling points


m = 4000;%The number of sampling nodes
g=1;%Derivatives of order at most 2
s1=logspace(-16,0,m).';%Sampling nodes
s2 = randi([0 g],m,1);%Random derivative for each j

S=blkdiag();%Jordan-like matrix using sampling nodes
for j=1:m
    S=blkdiag(S,diag(s1(j)*ones(s2(j)+1,1))+diag(ones(1,s2(j)),1));
end

y1=zeros(m,g+1);
y1(1:m,1) = s1.*sqrt(s1);%Function values at sampling nodes
y1(1:m,2) = 3/2*sqrt(s1);%First derivative at sampling nodes

yy=[];%Constructing the exact solution as the structure of the starting vector v2
yy1=[];
for j=1:m
    if s2(j)==0
        yy=[yy;y1(j,1)];
    else
        yy1=[yy1;y1(j,2)];
        yy=[yy;y1(j,1)];
    end
end


M1=2000;%The number of nodes
z=logspace(-16,0,M1).';%Nodes
w=ones(1,M1);%Weights
s = randi([0 g],M1,1);%Random derivative for each node



f=zeros(M1,g+1);
f(1:M1,1) = z.*sqrt(z);%Function values at nodes
f(1:M1,2) = 3/2*sqrt(z);%First derivative at nodes


ff1=[];%Constructing the right hand side vector the same as the structure of the starting vector v2
for j=1:M1
    ff=zeros(s(j)+1,1);
    for i=1:s(j)+1
        ff(i,1)=f(j,s(j)+1-(i-1));
        ff1=[ff1;ff(i,1)];
    end
end

v=[];%The starting vector v2 and the vector v
v2 = [];
for j=1:M1
    ww=zeros(s(j)+1,1);
    vv=zeros(s(j)+1,1);
    for i=1:s(j)+1
        ww(i,1)=w(1,j);
        if i==s(j)+1
            vv(i,1)=w(1,j);
        else
            vv(i,1)=0;
        end
        v=[v;ww(i,1)];
        v2=[v2;vv(i,1)];
    end
end

J=blkdiag();%Jordan-like matrix for the nodes
for j=1:M1
    J=blkdiag(J,diag(z(j)*ones(s(j)+1,1))+diag(ones(1,s(j)),1));
end

sum=0;%The Sobolev orthonormal rational function of degree zero
for k=1:M1
    sum=sum+w(k)^2;
end
r0=1/sqrt(sum);



max_degree = 100;%The maximum degree of the approximation
err1 = zeros(max_degree,1);
err2 = zeros(max_degree,1);


err21 = zeros(max_degree,1);
err22 = zeros(max_degree,1);


err31 = zeros(max_degree,1);
err32 = zeros(max_degree,1);

for n = 1:length(err1)
    C = 2;%Tapered lightening poles
    sigma = sqrt(2)*pi;
    pp = -C*exp(-sigma*(sqrt(n)-sqrt(1:n)));
    t=pp;%The vector of poles



    [y, H, K] = ratfitJ(J, t, ff1, n, v, v2);%Arnoldi-based approach
    r = ratvalJ(y, H, K, S, s2, r0);


    d2 = ratfitf1(z,s,t,ff1,n);%The approach of not employing Arnoldi
    y4 = ratvalf1(d2,s1,s2,t,n);


    d1 = ratfitf(z,s,t,ff1,n);%The method introduced in [15]
    y3 = ratvalf(d1,s1,s2,t,n);

    xx=[];%Separating the approximation of the function and its first derivative
    xx1=[];
    s3=[-1;s2];
    sum3=0;
    for j=2:m+1
        if j==2
            sum3=sum3+s3(j-1);
        else
            sum3=sum3+s3(j-1)+1;
        end
        if s3(j)==0
            xx=[xx;r(sum3+1+1,1)];
        else
            xx1=[xx1;r(sum3+1+1,1)];
            xx=[xx;r(sum3+1+2,1)];
        end
    end


    xx21=[];%Separating the approximation of the function and its derivatives for the approach without Arnoldi
    xx22=[];
    s4=[-1;s2];
    sum4=0;
    for j=2:m+1
        if j==2
            sum4=sum4+s4(j-1);
        else
            sum4=sum4+s4(j-1)+1;
        end
        if s4(j)==0
            xx21=[xx21;y3(sum4+1+1,1)];
        else
            xx22=[xx22;y3(sum4+1+1,1)];
            xx21=[xx21;y3(sum4+1+2,1)];
        end
    end

    xx31=[];%Separating the approximation of the function and its derivatives for the method of [15]
    xx32=[];
    s5=[-1;s2];
    sum5=0;
    for j=2:m+1
        if j==2
            sum5=sum5+s5(j-1);
        else
            sum5=sum5+s5(j-1)+1;
        end
        if s5(j)==0
            xx31=[xx31;y4(sum5+1+1,1)];
        else
            xx32=[xx32;y4(sum5+1+1,1)];
            xx31=[xx31;y4(sum5+1+2,1)];
        end
    end



    err1(n) = max(abs(xx - yy));%Errors for Arnoldi
    err2(n) = max(abs(xx1 - yy1));


    err21(n) = max(abs(xx31 - yy));%Errors for the approach without Arnoldi
    err22(n) = max(abs(xx32 - yy1));



    err31(n) = max(abs(xx21 - yy));%Errors for the approach introduced in [15]
    err32(n) = max(abs(xx22 - yy1));


end
%The figure of the error of r
plot(sqrt(1:length(err1)), log10(err1),"b.",sqrt(1:length(err21)), log10(err21),"ro",sqrt(1:length(err31)), log10(err31),"g*")
legend('Confluent Cauchy with Arnoldi','Confluent Cauchy without Arnoldi 1','Confluent Cauchy without Arnoldi 2','Location','southwest')

xlabel('$\sqrt{n}$','Interpreter','Latex');
ylabel('$Log_{10}(||f-r||_{\infty})$','Interpreter','Latex');

saveas(gcf, 'exm41.png');

%The figure of the error of r'
% plot(sqrt(1:length(err2)), log10(err2),"b.",sqrt(1:length(err22)), log10(err22),"ro",sqrt(1:length(err32)), log10(err32),"g*")
% legend('Confluent Cauchy with Arnoldi','Confluent Cauchy without Arnoldi 1','Confluent Cauchy without Arnoldi 2','Location','southwest')
%
% xlabel('$\sqrt{n}$','Interpreter','Latex');
% ylabel('$Log_{10}(||f^{\prime}-r^{\prime}||_{\infty}))$','Interpreter','Latex');
%
% saveas(gcf, 'exm42.png');

%Arnoldi-based approach
function [y,H,K] = ratfitJ(J,Pole,f,n,w,v)%The function related to Algorithm 7
m = length(v);
Q = v/norm(v);
H = zeros(n+1,n);
for k = 1:n
    if k == 1
        q = ((J-Pole(1)*eye(m))\eye(m))*Q(:,1);
    else
        q = ((J-Pole(k)*eye(m))\(J-Pole(k-1)*eye(m)))*Q(:,k);
    end
    for j = 1:k
        H(j,k) = Q(:,j)'*(q);
        q = q - H(j,k)*Q(:,j);
    end
    H(k+1,k) = norm(q);
    Q = [Q q/H(k+1,k)];
end
E = eye(n+1,n); E(1,1) = 0;
K = H - E;
H = H*diag(Pole) - eye(n+1,n)*diag([-1,Pole(1:n-1)]);
y = Q\(w.*f);
end
%Approximation
function r = ratvalJ(y,H,K,X,S,r0)%The function related to Algorithm 8
tau = length(S);
M = size(X,1);
n = size(H,2);
U = [];
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
r = U*y;
end

%The approach of not employing Arnoldi
function d2 = ratfitf1(z,s,t,ff1,n)
M1=length(z);
A1=zeros(M1,n+1);
A1(:,1)=ones(M1,1);
for i=1:n
    A1(:,i+1)=1./(z-t(i));
end

A2=zeros(M1,n+1);
A2(:,1)=zeros(M1,1);
for i=1:n
    A2(:,i+1)=-1./(z-t(i)).^2;
end
A=[];%Constructing the confluent Cauchy-like matrix
for j=1:M1
    if s(j)==0
        A=[A;A1(j,:)];
    else
        A=[A;A2(j,:);A1(j,:)];
    end
end
d2=A\ff1;%Solving the system including a confluent Cauchy-like matrix directly
end
%Approximation
function y4 = ratvalf1(d2,s1,s2,t,n)
m=length(s1);
B1=zeros(m,n+1);
B1(:,1)=ones(m,1);
for i=1:n
    B1(:,i+1)=1./(s1-t(i));
end

B2=zeros(m,n+1);
B2(:,1)=zeros(m,1);
for i=1:n
    B2(:,i+1)=-1./(s1-t(i)).^2;
end
B=[];%Constructing the confluent Cauchy-like matrix for sampling points
for j=1:m
    if s2(j)==0
        B=[B;B1(j,:)];
    else
        B=[B;B2(j,:);B1(j,:)];
    end
end
y4=B*d2;%The approximation by solving the system including a confluent Cauchy-like matrix directly
end

%The method introduced in [15]
function d1 = ratfitf(z,s,t,ff1,n)
M1=length(z);
A1=zeros(M1,n+1);
A1(:,1)=ones(M1,1);
for i=1:n
    A1(:,i+1)=t(i)./(z-t(i));
end

A2=zeros(M1,n+1);
A2(:,1)=zeros(M1,1);
for i=1:n
    A2(:,i+1)=-t(i)./(z-t(i)).^2;
end
A=[];%Constructing the confluent Cauchy-like matrix with poles in the numerators (method introduced in [15])
for j=1:M1
    if s(j)==0
        A=[A;A1(j,:)];
    else
        A=[A;A2(j,:);A1(j,:)];
    end
end

d1=A\ff1;
end

%Approximation
function y3 = ratvalf(d1,s1,s2,t,n)
m=length(s1);
B1=zeros(m,n+1);
B1(:,1)=ones(m,1);
for i=1:n
    B1(:,i+1)=t(i)./(s1-t(i));
end

B2=zeros(m,n+1);
B2(:,1)=zeros(m,1);
for i=1:n
    B2(:,i+1)=-t(i)./(s1-t(i)).^2;
end
B=[];%Constructing the Cauchy-like matrix with poles in the numerators for sampling nodes (method introduced in [15])
for j=1:m
    if s2(j)==0
        B=[B;B1(j,:)];
    else
        B=[B;B2(j,:);B1(j,:)];
    end
end
y3=B*d1;%The approximation by solving the system including a confluent Cauchy-like matrix directly
end
