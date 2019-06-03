function [c,f,s] = pde_1_function(x,t,u,DuDx)
c = 1;
f = 0.47*DuDx;
s = - 1.55*u*DuDx;