function [u,x,y,t]=ADI_method(T,K,N,M)
% solves the 2D heat equation u_{t}=K u_{xx}+ f(x,y,t) in the unit square 
% and on the interval [0,T] with zero boundary condition by ADI method.
%K is the (heat conduction) coefficient in the heat equation.
%N is the number of subintervals in [0,1].
%M is the number of time steps.
h=1/N;
tau=T/M;
gamma=tau*K/h^2;
u=zeros(N+1,N+1,M+1);
u1=zeros(N+1,N+1);
alpha=zeros(N,1);
beta=zeros(N,1);
% grid points
x=(0:N)*h;
y=x;
t=(0:M)*tau;

%Split the problem into the homogeneous PDE/non-homogeneous PDE parts as stated
%in the document.

% want to find the solution v where:
% v(x,y,t) = u(x,y,t) + w(x,y,t)

% u - zero BC, non-homogeneous PDE
% w - non-zero BC, homogeneous PDE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Approxiamte u, using the ADI method with Double Sweep   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initial condition
for k=2:N
    for j=2:N
        u(k,j,1)=initial_cond(x,y);
    end
end
% ADI method using the double-sweep method
for k=1:N-1
    alpha(k+1)=gamma/(2*(1+gamma)-gamma*alpha(k));
end
for n=2:M+1
    for j=2:N
        for k=2:N
            F=-0.5*gamma*(u(k,j-1,n-1)+u(k,j+1,n-1))-(1-gamma)*u(k,j,n-1)...
                -0.25*tau*(ff(x(k),y(j),t(n))+ff(x(k),y(j),t(n-1)));
            beta(k)=(beta(k-1)*gamma-2*F)/(2*(1+gamma)-gamma*alpha(k-1));
        end
        for k=N:(-1):2
            u1(k,j)=alpha(k)*u1(k+1,j)+beta(k);
        end
    end
    for k=2:N
        for j=2:N
            F=-0.5*gamma*(u1(k-1,j)+u1(k+1,j))-(1-gamma)*u1(k,j)...
                -0.25*tau*(ff(x(k),y(j),t(n))+ff(x(k),y(j),t(n-1)));
            beta(j)=(beta(j-1)*gamma-2*F)/(2*(1+gamma)-gamma*alpha(j-1));
        end
        for j=N:(-1):2
            u(k,j,n)=alpha(j)*u(k,j+1,n)+beta(j);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   Compute w at the grid points for each t  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose w(x,y,t) = x*y^2 -2y -2x + y^x^2, which satisfies the BCs.

w = zeros(N+1,N+1);
for k=1:N+1
    for j=1:N+1
        w(k,j) = x(k)^2*y(j)^2 - 2*x(k) - 2*y(j);
    end
end

% Add our solutions togethet to get our final approximation.

for n=1:M+1
    u(:,:,n) = u(:,:,n) + w;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function f(x,y,t) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=ff(x,y,t)
f=24*x^2*(x^2 + y^2) + 2*0.52*(x^2 + y^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function g(x,y) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g=initial_cond(x,y)
g=0;
        