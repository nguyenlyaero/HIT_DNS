function [fx,fy,fz] = project_divergence_free(f1,f2,f3)
[nx,ny,nz] = size(f1);
kx = fftshift(-nx/2:1:nx/2-1); ky = fftshift(-ny/2:1:ny/2-1); kz = fftshift(-nz/2:1:nz/2-1);
[k1,k2,k3] = ndgrid(kx,ky,kz);

onebyksq = 1./(k1.*k1 + k2.*k2 + k3.*k3 + 1E-16);
onebyksq(1,1,1) = 0;

f1hat = fftn(f1);
f2hat = fftn(f2);
f3hat = fftn(f3);

fx = real(ifftn( (1 - k1.*k1.*onebyksq).*f1hat - (k1.*k2.*onebyksq).*f2hat - (k1.*k3.*onebyksq).*f3hat));
fy = real(ifftn(-(k2.*k1.*onebyksq).*f1hat + (1 - k2.*k2.*onebyksq).*f2hat - (k2.*k3.*onebyksq).*f3hat));
fz = real(ifftn(-(k3.*k1.*onebyksq).*f1hat - (k3.*k2.*onebyksq).*f2hat + (1 - k3.*k3.*onebyksq).*f3hat));

