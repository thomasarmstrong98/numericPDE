function [pl,ql,pr,qr] = pde_1_bc(xl,ul,xr,ur,t)
pl = ul - 0.5;
ql = 0;
pr = 0;
qr = 1;