clear; clc;
addpath('/scratch/06005/nguyenly/HIT_DNS/PadeOps_output');
addpath('../MATLAB');

Run=5;
N=512;
Re=114.85;
N0=60;
Nsample=39;
% Read from Last Steady Restart File
u0 = read_fortran_box(['RESTART_Run0' num2str(Run, '%d') '_u.00' num2str(N0,'%02d') '00'], N, N, N, 'double');
v0 = read_fortran_box(['RESTART_Run0' num2str(Run, '%d') '_v.00' num2str(N0,'%02d') '00'], N, N, N, 'double');
w0 = read_fortran_box(['RESTART_Run0' num2str(Run, '%d') '_w.00' num2str(N0,'%02d') '00'], N, N, N, 'double');
fid =fopen(['RESTART_Run0' num2str(Run, '%d') '_info.00' num2str(N0,'%02d') '00']);
t=fscanf(fid,'%f'); t0=t(1);

% Enstrophy Spectrum
omega_x=ddy_hit(w0)-ddz_hit(v0);
omega_y=ddz_hit(u0)-ddx_hit(w0);
omega_z=ddx_hit(v0)-ddy_hit(u0);

[Omega_x,~]=get_energy_spectrum(omega_x,500);
[Omega_y,~]=get_energy_spectrum(omega_y,500);
[Omega_z,kline]=get_energy_spectrum(omega_z,500);
kline=kline.';

Omega0=(Omega_x+Omega_y+Omega_z).^(1/2);

% Energy Spectrum
[Eu,~] = get_energy_spectrum(u0,500);
[Ev,~] = get_energy_spectrum(v0,500);
[Ew,~] = get_energy_spectrum(w0,500);

E0=1/2*(Eu+Ev+Ew);

% Mean TKE
E=(u0.*u0 + v0.*v0 + w0.*w0)/2;
Emean0=mean(E(:));

% Energy Dissipation Rate
dudx=ddx_hit(u0); dudy=ddy_hit(u0); dudz=ddz_hit(u0);
dvdx=ddx_hit(v0); dvdy=ddy_hit(v0); dvdz=ddz_hit(v0);
dwdx=ddx_hit(w0); dwdy=ddy_hit(w0); dwdz=ddz_hit(w0);

epsilon0=dudx.*dudx+dudy.*dudy+dudz.*dudz...
    +dvdx.*dvdx+dvdy.*dvdy+dvdz.*dvdz...
    +dwdx.*dwdx+dwdy.*dwdy+dwdz.*dwdz;
    
epsilon0=1/Re*mean(epsilon0(:));

EMat=zeros(Nsample,500);
OmegaMat=zeros(Nsample,500);
tvec=zeros(Nsample,1);
Evec=zeros(Nsample,1);
epsilonvec=zeros(Nsample,1);

for n=1:Nsample
u = read_fortran_box(['Run0' num2str(Run+3, '%d') '_uVel_t00' num2str(n+N0,'%02d') '00.out'], N, N, N, 'double');
v = read_fortran_box(['Run0' num2str(Run+3, '%d') '_vVel_t00' num2str(n+N0,'%02d') '00.out'], N, N, N, 'double');
w = read_fortran_box(['Run0' num2str(Run+3, '%d') '_wVel_t00' num2str(n+N0,'%02d') '00.out'], N, N, N, 'double');
fid =fopen(['Run0' num2str(Run+3, '%d') '_info_t00' num2str(n+N0,'%02d') '00.out']);
t=fscanf(fid,'%f'); t=t(1);
tvec(n)=t-t0;

omega_x=ddy_hit(w)-ddz_hit(v);
omega_y=ddz_hit(u)-ddx_hit(w);
omega_z=ddx_hit(v)-ddy_hit(u);

[Omega_x,~]=get_energy_spectrum(omega_x,500);
[Omega_y,~]=get_energy_spectrum(omega_y,500);
[Omega_z,~]=get_energy_spectrum(omega_z,500);

OmegaMat(n,:)=(Omega_x+Omega_y+Omega_z).^(1/2);

[Eu,~] = get_energy_spectrum(u,500);
[Ev,~] = get_energy_spectrum(v,500);
[Ew,~] = get_energy_spectrum(w,500);

EMat(n,:)=1/2*(Eu+Ev+Ew);

% Kinetic Energy
E=(u.*u + v.*v + w.*w)/2;
E=mean(E(:));
Evec(n)=E;

% Energy dissipation rate
dudx=ddx_hit(u); dudy=ddy_hit(u); dudz=ddz_hit(u);
dvdx=ddx_hit(v); dvdy=ddy_hit(v); dvdz=ddz_hit(v);
dwdx=ddx_hit(w); dwdy=ddy_hit(w); dwdz=ddz_hit(w);

epsilon=dudx.*dudx+dudy.*dudy+dudz.*dudz...
    +dvdx.*dvdx+dvdy.*dvdy+dvdz.*dvdz...
    +dwdx.*dwdx+dwdy.*dwdy+dwdz.*dwdz;
    
epsilon=1/Re*mean(epsilon(:));
epsilonvec(n)=epsilon;

% Print
fprintf('%d \n', n);


end

% Plot Avg Energy Spectrum
figure;
title('Energy Spectra');
loglog(kline, E0);
hold on;
for i=1:Nsample
    loglog(kline, EMat(i,:));
end
ylim([1e-18 1e2]);
xlim([0.5, inf]);
legendCell=cellstr(num2str(tvec', 't=%f'));
legendCell={'t=0', legendCell};
legend(legendCell);
hold off;
saveas(gcf,'ESpectraDecay.fig');
close;
        
% Plot Avg Enstrophy Spectrum
figure;
title('Ensthophy Spectra');
loglog(kline, Omega0);
hold on;
for i=1:Nsample
    loglog(kline, OmegaMat(i,:));
end
ylim([1e-18 1e2]);
xlim([0.5, inf]);
legendCell=cellstr(num2str(tvec', 't=%f'));
legendCell={'t=0', legendCell};
legend(legendCell);
hold off;
saveas(gcf,'OmegaSpectraDecay.fig');
close;

% Plot TKE vs scaled time
tvec=[0; tvec];
Evec=[Eman0; Evec];
epsilonvec=[epsilon0; epsilonvec];
figure;
title('TKE');
plot(tvec, Evec);
xlabel('t');
ylabel('TKE');
saveas(gcf,'TKEDecay.fig');
close;

save('decay.mat', 'EMat', 'kline', 'OmegaMat', 'tvec', 'Evec', 'epsilonvec');
