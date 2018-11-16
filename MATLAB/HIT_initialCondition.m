clear, clc, close all

seed = 12345;
nx = 512; ny = 512; nz = 512;
[u,v,w] = generate_isotropic_turbulence(nx,ny,nz,seed);

dat = zeros(nx*ny*nz,3);
dat(:,1) = real(u(:));
dat(:,2) = real(v(:));
dat(:,3) = real(w(:));

save('PadeOps_HIT.dat','dat','-ascii')