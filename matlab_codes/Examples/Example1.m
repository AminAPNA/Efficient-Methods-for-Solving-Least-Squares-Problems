%This code is done for Example 1, in which a Sobolev polynomial least squares problem is considered
%It Generates three plots related to Example 1 (Figure 1)
%INPUT:
%   J = MxM matrix (We consider M nodes)
%   v = the column vector of order M including weights and factors compensating for the repeated differentiation which introduces factorial terms
%   v2 = column vector of size M
%   ff1 = the column vector including the data points
%   X = the Jordan-like matrix constructed by the sampling nodes s1
%   s2 = The vector of order m including the sampling derivatives \tilde{s}_{j}
%   p0 = the Sobolev orthonormal polynomial of degree 0 directly computed through the inner product


%First we compute the matrix Q, H and also vector y through the function polyfitJ
%Then by taking a set of sampling nodes, we evaluate the approximations by polyvalJ.
%At last we compute the errors at sampling points


m = 1000;%The number of sampling nodes
g=2;%Derivatives of order at most 2
[s1]=jacpts(m,-.5,-.5);%Sampling nodes
s2 = randi([0 g],m,1);%Sampling derivatives for each j

X=blkdiag();%Jordan-like matrix using sampling nodes
for j=1:m
    if s2(j)==0 || s2(j)==1
        X=blkdiag(X,diag(s1(j)*ones(s2(j)+1,1))+diag(ones(1,s2(j)),1));
    else
        X=blkdiag(X,diag(s1(j)*ones(s2(j)+1,1))+diag([2 1],1));
    end
end

y1=zeros(m,g+1);
y1(1:m,1) = 1./(1+25*s1.^2);%Function values at sampling nodes
y1(1:m,2) = -50*s1./((1+25*s1.^2).^2); %First derivative at sampling nodes
y1(1:m,3) = ((-50.*(1+25*s1.^2))+50*100*s1.^2)./((1+25*s1.^2).^3);%Second derivative at sampling nodes


yy=[];%Constructing the exact solution and derivatives as the structure of the starting vector of the Arnoldi iteration v2
yy1=[];
yy2=[];
for j=1:m
    if s2(j)==0
        yy=[yy;y1(j,1)];
    elseif s2(j)==1
        yy1=[yy1;y1(j,2)];
        yy=[yy;y1(j,1)];
    else
        yy2=[yy2;y1(j,3)];
        yy1=[yy1;y1(j,2)];
        yy=[yy;y1(j,1)];
    end
end

max_degree = 200;%The maximum degree of the approximation n
err1 = zeros(max_degree,1);%Error vectors for the approach with Arnoldi
err2 = zeros(max_degree,1);
err3 = zeros(max_degree,1);


err11 = zeros(max_degree,1);%Error vectors for the approach without Arnoldi
err21 = zeros(max_degree,1);
err31 = zeros(max_degree,1);

for n = 1:length(err1)%The loop for different n
    M=2*n+1;%The number of nodes
    [z,w]=jacpts(M,-.5,-.5);%The pair of node-weight
    s = randi([0 g],M,1);%Random derivative for each node


    J=blkdiag();%Jordan-like matrix for the nodes
    for j=1:M
        if s(j)==0 || s(j)==1
            J=blkdiag(J,diag(z(j)*ones(s(j)+1,1))+diag(ones(1,s(j)),1));
        else
            J=blkdiag(J,diag(z(j)*ones(s(j)+1,1))+diag([2 1],1));
        end
    end


    f=zeros(M,g+1);
    f(1:M,1) = 1./(1+25*z.^2);%Function values at nodes
    f(1:M,2) = -50*z./((1+25*z.^2).^2); %First derivative at nodes
    f(1:M,3) = ((-50.*(1+25*z.^2))+50*100*z.^2)./((1+25*z.^2).^3);%Second derivative at nodes



    ff1=[];%Constructing the right hand side vector the same as the structure of the starting vector v2
    for j=1:M
        ff=zeros(s(j)+1,1);
        for i=1:s(j)+1
            ff(i,1)=f(j,s(j)+1-(i-1));
            ff1=[ff1;ff(i,1)];
        end
    end

    v=[];%The starting vector v2 and the vector v
    v2 = [];
    for j=1:M
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

    sum=0;%Polynomial of degree zero
    for k=1:M
        sum=sum+w(k)^2;
    end
    p0=1/sqrt(sum);


    [y, H] = polyfitJ(J, ff1, n, v, v2);%Arnoldi-based approach
    p = polyvalJ(y, H, X, s2, p0);


    d1 = polyfitf(z,s,v,ff1,n);%Whithout Arnoldi
    y1 = polyvalf(d1,s1,s2,n);


    xx=[];%Separating the approximation of the function and its derivatives
    xx1=[];
    xx2=[];
    s3=[-1;s2];
    sum3=0;
    for j=2:m+1
        if j==2
            sum3=sum3+s3(j-1);
        else
            sum3=sum3+s3(j-1)+1;
        end
        if s3(j)==0
            xx=[xx;p(sum3+1+1,1)];
        elseif s3(j)==1
            xx1=[xx1;p(sum3+1+1,1)];
            xx=[xx;p(sum3+1+2,1)];
        else
            xx2=[xx2;p(sum3+1+1,1)];
            xx1=[xx1;p(sum3+1+2,1)];
            xx=[xx;p(sum3+1+3,1)];
        end
    end

    xx21=[];%Separating the approximation of the function and its derivatives for without Arnoldi
    xx22=[];
    xx23=[];
    s4=[-1;s2];
    sum4=0;
    for j=2:m+1
        if j==2
            sum4=sum4+s4(j-1);
        else
            sum4=sum4+s4(j-1)+1;
        end
        if s4(j)==0
            xx21=[xx21;y1(sum4+1+1,1)];
        elseif s4(j)==1
            xx22=[xx22;y1(sum4+1+1,1)];
            xx21=[xx21;y1(sum4+1+2,1)];
        else
            xx23=[xx23;y1(sum4+1+1,1)];
            xx22=[xx22;y1(sum4+1+2,1)];
            xx21=[xx21;y1(sum4+1+3,1)];
        end
    end


    err1(n) = max(abs(xx - yy));%Errors for Arnoldi
    err2(n) = max(abs(xx1 - yy1));
    err3(n) = max(abs(xx2 - yy2));

    err11(n) = max(abs(xx21 - yy));%Errors for the approach without Arnoldi
    err21(n) = max(abs(xx22 - yy1));
    err31(n) = max(abs(xx23 - yy2));
end
%The plot for the errors of p
plot(1:length(err1), log10(err1(1:200)),"b.",1:length(err11), log10(err11(1:200)),"ro")%Figure of Function approximation errores
legend('Confluent Vandermonde with Arnoldi','Confluent Vandermonde','Location','southwest')

xlabel('$n$','Interpreter','Latex');
ylabel('$Log_{10}(||f-p||_{\infty})$','Interpreter','Latex');

saveas(gcf, 'exm21.png');

%The plot for the errors of p'
%   plot(1:length(err2), log10(err2(1:200)),"b.",1:length(err21), log10(err11(1:200)),"ro")%Figure of the first derivative approximation errores
%  legend('Confluent Vandermonde with Arnoldi','Confluent Vandermonde','Location','southwest')
% 
% xlabel('$n$','Interpreter','Latex');
% ylabel('$Log_{10}(||f^{\prime}-p^{\prime}||_{\infty})$','Interpreter','Latex');
% 
% saveas(gcf, 'exm22.png');

% The plot for the errors of p"
   plot(1:length(err3), log10(err3(1:200)),"b.",1:length(err31), log10(err11(1:200)),"ro")%Figure of the second derivative approximation errores
 legend('Confluent Vandermonde with Arnoldi','Confluent Vandermonde','Location','southwest')
xlabel('$n$','Interpreter','Latex');
ylabel('$Log_{10}(||f^{\prime \prime}-p^{\prime \prime}||_{\infty})$','Interpreter','Latex');
saveas(gcf, 'exm23.png');


%Arnoldi-based approach, the function related to Algorithm 3 of the paper
function [y,H] = polyfitJ(J,f,n,w,v)
Q = v/norm(v);
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
y = Q\(w.*f);
end
%Approximation, the function related to Algorithm 4 of the paper
function p = polyvalJ(y,H,X,S,p0)
tau = length(S);
n = size(H,2);
U = [];
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
p = U*y;
end

%Without Arnoldi approach
function d1 = polyfitf(z,s,v,ff1,n)%This function is provided to compute the vector of unknown directly by solving a confluent Vandermonde system
M=length(z);
A1=z.^(0:n);
A2=zeros(M,n);
for i=1:M
    for j=1:n
        A2(i,j)=(j)*(z(i)^(j-1));
    end
end
A2=[zeros(M,1) A2(:,1:n)];

A3=zeros(M,n-1);
for i=1:M
    for j=1:n
        A3(i,j)=(j)*(j+1)*(z(i)^(j-1));
    end
end
A3=[zeros(M,1) zeros(M,1) A3(:,1:n-1)];

A=[];%Constructing the confluent Vandermonde matrix which includes the monomials and its derivatives at nodes
for j=1:M
    if s(j)==0
        A=[A;A1(j,:)];
    elseif s(j)==1
        A=[A;A2(j,:);A1(j,:)];
    else
        A=[A;A3(j,:);A2(j,:);A1(j,:)];
    end
end
A=diag(v)*A;
d1=A\(v.*ff1);
end
%Approximation
function y1 = polyvalf(d1,s1,s2,n)%This function is provided to compute the approximations by solving a confluent Vandermonde system directly
b=length(s1);
B1=s1.^(0:n);
B2=zeros(b,n);
for i=1:b
    for j=1:n
        B2(i,j)=(j)*(s1(i)^(j-1));
    end
end
B2=[zeros(b,1) B2(:,1:n)];

B3=zeros(b,n-1);
for i=1:b
    for j=1:n
        B3(i,j)=(j)*(j+1)*(s1(i)^(j-1));
    end
end
B3=[zeros(b,1) zeros(b,1) B3(:,1:n-1)];

B=[];%Constructing the confluent Vandermonde matrix which includes the monomials and its derivatives at sampling nodes
for j=1:b
    if s2(j)==0
        B=[B;B1(j,:)];
    elseif s2(j)==1
        B=[B;B2(j,:);B1(j,:)];
    else
        B=[B;B3(j,:);B2(j,:);B1(j,:)];
    end
end
y1=B*d1;%Solving the system including the confluent Vandermonde matrix directly (Whithout Arnoldi)
end