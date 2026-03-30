%This code is done for Example 2, in which a rational least squares problem is considered
%It Generates Figure 2 related to Example 2
%INPUT:
%   z = The vector of order M1, including the nodes
%   t = The column vector including n poles
%   w = column vector of size M1 including weights
%   f = the column vector including the data values
%   s = the column vector s of order m including the sampling nodes
%   r0 = the orthonormal rational function of degree 0 directly computed through the inner product


%First we compute the matrix Q, H, K and also vector y through the function ratfitZ
%Then by taking a set of sampling nodes, we evaluate the approximations by ratvalZ.
%At last we compute the errors at sampling points


m = 2000;%Sampling nodes
s=logspace(-16,0,m).';
s = [s; -s]; m = 2*m;
y1=abs(s);%Exact function values at sampling nodes

M1=1000;%Nodes and weights
z=logspace(-16,0,M1).';%Exponentially clustered nodes
z = [z; -z]; M1 = M1*2;
w=ones(1,M1);
w=w';

f=abs(z);%Function values at nodes

sum=0;%The first rational function
for k=1:M1
    sum=sum+w(k)^2;
end
r0=1/sqrt(sum);


max_degree = 100;%Max number of n
err = zeros(max_degree,1);
err1 = zeros(max_degree,1);
err2 = zeros(max_degree,1);


for n = 1:length(err)
    C = 2; %Tapered lightening poles
    sigma = sqrt(2)*pi;
    pp = 1i*sqrt(C*exp(-sigma*(sqrt(n)-sqrt(0:n))));
    pp = [pp, -pp];
    t=pp;%The vector of poles

    [y, H, K] = ratfitZ(z, t, f, 2*(n+1), w);%Rational Arnoldi
    r = ratvalZ(y, H, K, s, r0);

    c = ratfit1(z, f, t, 2*(n+1));%The approach of not employing Arnoldi
    y2 = ratval1(c, s, t);

    c1 = ratfit2(z, f, t, 2*(n+1));%The method introduced in [15]
    y3 = ratval2(c1, s, t);


    err(n) = max(abs(r - y1));%Errors
    err1(n) = max(abs(y2 - y1));
    err2(n) = max(abs(y3 - y1));
end
%The figure of errors
plot(sqrt(1:length(err)), log10(err),"b.",sqrt(1:length(err1)), log10(err1),"ro",sqrt(1:length(err2)), log10(err2),"g*");
legend('Arnoldi-based approach','The approach of not employing Arnoldi','The method introduced in [15]','Location','southwest')
xlabel('$\sqrt{n}$','Interpreter','Latex');
ylabel('$Log_{10}(||f-r||_{\infty})$','Interpreter','Latex');
saveas(gcf, 'exm31.png');
export_fig test2.png

%Arnoldi-based approach
function [y,H,K] = ratfitZ(z,Pole,f,n,v)%The function related to Algorithm 5
m = length(z);
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
    Q = [Q q/H(k+1,k)];
end
E = eye(n+1,n); E(1,1) = 0;
K = H - E;
H = H*diag(Pole) - eye(n+1,n)*diag([-1,Pole(1:n-1)]);
y = Q\(v.*f);
end

function r = ratvalZ(y,H,K,x,r0)%The function related to Algorithm 6
M = length(x);
n = size(H,2);
U = zeros(M,1);
U(:,1) = r0;
for k = 1:n
    u = 0;
    for j = 1:k
        u = u + H(j,k)*U(:,j) - K(j,k)*x.*U(:,j);
    end
    U = [U u./(K(k+1,k)*x - H(k+1,k))];
end
r = U*y;
end

%The approach of not employing Arnoldi
function c = ratfit1(z,f,t,n)
M1=length(z)-1;
A=zeros(M1+1,n+1);%Constructing the Cauchy-like matrix
for i=1:n
    A(:,1)=ones(M1+1,1);
    A(:,i+1)=1./(z-t(i));
end
c=A\f;%Solving the system including a Cauchy-like matrix directly
end

function y2 = ratval1(c,s,t)
m=length(s);
n=length(c)-1;
B=zeros(m,n+1);%Constructing the Cauchy-like matrix for sampling points
for i=1:n
    B(:,1)=ones(m,1);
    B(:,i+1)=1./(s'-t(i));
end
y2=B*c;%The approximation by solving the system including a Cauchy-like matrix directly
end

%The method introduced in [15]
function c1 = ratfit2(z,f,t,n)
M1=length(z)-1;
A=zeros(M1+1,n+1);%Constructing the Cauchy-like matrix with poles in the numerators (method introduced in [15])
for i=1:n
    A(:,1)=ones(M1+1,1);
    A(:,i+1)=t(i)./(z-t(i));
end
c1=A\f;
end

function y3 = ratval2(c,s,t)
m=length(s);
n=length(c)-1;
B=zeros(m,n+1);%Constructing the Cauchy-like matrix with poles in the numerators for sampling nodes (method introduced in [15])
for i=1:n
    B(:,1)=ones(m,1);
    B(:,i+1)=t(i)./(s'-t(i));
end
y3=B*c;%The approximation by solving the system including a Cauchy-like matrix directly
end