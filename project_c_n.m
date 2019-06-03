function u = project_c_n(T,M,x,u0)
% Solves the heat equation by Cranc-Nicolson method;
% T - specifies the last time point we want to evaluate
% M - M+1, the number of discrete time points we consider for time
% x - (pre-calculated) array of grid points in x, 
% u0 - (pre-calculated) initial condition at grid points
N=length(x)-1;
if (N+1)~=length(u0)
    error('Dimensions of x and u0 must agree')
end
% Since this problem has non-zero BCs, we will split into 2 functions
% u(x,t) = v(x,t) + g(x,t)
% v(x,t) is the solution with zero boundary conditions and
% g(x,t) is the solution which contains the boundary condition 
% g(x,t) = 0.9*x

% The inital condition u0 passed to this function is that of u(x,t), so we must modify
% it to make a v0. This is done by:
% v(x, 0) = u(x,0) - g(x,0) = 0.9*(sin(pi*x/2) - x)

g = @(x) 0.9*x;
v0 = u0 - g(x);
h=(x(N+1)-x(1))/N;
tau=T/M;
gamma=tau/h^2;
w=zeros(N+1,M+1);
% grid points in time
t=(0:M)*tau;
% initial condition
w(:,1)=v0;
% solving the tridiagonal system by the double-sweep method
alpha=zeros(1,N);
beta=zeros(1,N);

f = @(x) 25*x*(1-x^2);

for k=1:(N-1)
    alpha(k+1)=0.5*gamma/(1+0.5*gamma*(2-alpha(k)));
end
for m=1:M
    % finding the alpha and betas
    for k=1:(N-1)
        F=-0.5*gamma*(w(k,m)+w(k+2,m))-(1-gamma)*w(k+1,m)-tau*(f(x(k+1)));
        beta(k+1)=(beta(k)*0.5*gamma-F)/(1+0.5*gamma*(2-alpha(k)));
    end
    % back propagation to find w
    for k=(N+1):(-1):3
        w(k-1,m+1)=alpha(k-1)*w(k,m+1)+beta(k-1);
    end
end

% Since w is the solution v(x,t) which has 0 boundary conditions, we must 
% add the function g(x,t) which contained the BC's to obtain u(x,t).
% u(x,t) = v(x,t) + g(x,t)
% g(x,t) = 0.9*x

g_x = g(x)';
g_mat = repmat(g_x, 1, M + 1);
u = w + g_mat;
% surface plot of u(x,t)
u = u';
surf(x,t,u);
xlabel('x');
ylabel('t');
zlabel('u');