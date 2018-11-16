clear; clc;
addpath('/scratch/06005/nguyenly/HIT_DNS/PadeOps_output');
addpath('/home1/06005/nguyenly/PadeOps/MATLAB');

% Read from Restart File
u0 = read_fortran_box('RESTART_Run03_u.000000', 256, 256, 256, 'double');
v0 = read_fortran_box('RESTART_Run03_v.000000', 256, 256, 256, 'double');
w0 = read_fortran_box('RESTART_Run03_w.000000', 256, 256, 256, 'double');

omega_x=ddy_hit(w0)-ddz_hit(v0);
omega_y=ddz_hit(u0)-ddx_hit(w0);
omega_z=ddx_hit(v0)-ddy_hit(u0);

[Omega_x,~]=get_energy_spectrum(omega_x,500);
[Omega_y,~]=get_energy_spectrum(omega_y,500);
[Omega_z,kline]=get_energy_spectrum(omega_z,500);
kline=kline.';

Omega0=(Omega_x+Omega_y+Omega_z).^(1/2);

Emean=zeros(78,500);
Omegamean=zeros(78,500);
tvec=zeros(78,1);
Evec=zeros(78,1);
epsilonvec=zeros(78,1);

parfor n=1:78
u = read_fortran_box(['Run03_uVel_t00' num2str(n,'%02d') '00.out'], 256, 256, 256, 'double');
v = read_fortran_box(['Run03_vVel_t00' num2str(n,'%02d') '00.out'], 256, 256, 256, 'double');
w = read_fortran_box(['Run03_wVel_t00' num2str(n,'%02d') '00.out'], 256, 256, 256, 'double');
fid =fopen(['Run03_info_t00' num2str(n,'%02d') '00.out']);
t=fscanf(fid,'%f'); t=t(1);
tvec(n)=t;

omega_x=ddy_hit(w)-ddz_hit(v);
omega_y=ddz_hit(u)-ddx_hit(w);
omega_z=ddx_hit(v)-ddy_hit(u);

[Omega_x,~]=get_energy_spectrum(omega_x,500);
[Omega_y,~]=get_energy_spectrum(omega_y,500);
[Omega_z,~]=get_energy_spectrum(omega_z,500);

Omegamean(n,:)=(Omega_x+Omega_y+Omega_z).^(1/2);

[Eu,~] = get_energy_spectrum(u,500);
[Ev,~] = get_energy_spectrum(v,500);
[Ew,~] = get_energy_spectrum(w,500);

Emean(n,:)=1/2*(Eu+Ev+Ew);

% if mod(n,10)==0
% 
% % Plot u
% figure, surface(real(u(:,:,45))','edgecolor','none'), daspect([1 1 1]);
% title(['t=' num2str(t)]);
% saveas(gcf,[num2str(n) '00u.pdf']);
% close;
% 
% % Plot energy spectrum
% figure;
% title(['t=' num2str(t)]);
% loglog(kline, 1/2*(Eu+Ev+Ew));
% hold on;
% loglog(kline, 1.8191*114.85/34*1.4528*(kline.^4)./((1 + kline.^2).^(17/6)).*exp(-5.2.*kline.*0.0245));
% ylim([1e-18 1e2]);
% xlim=([0.4 400]);
% hold off;
% saveas(gcf,[num2str(n) '00E.fig']);
% close;
% 
% % Plot enstrophy spectrum
% figure;
% title(['t=' num2str(t)]);
% loglog(klineO, Omegamean(n,:));
% hold on;
% loglog(klineO, Omega0);
% hold off;
% saveas(gcf,[num2str(n) '00Omega.fig']);
% close;
% 
% end

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
    
epsilon=1/114.85*mean(epsilon(:));
epsilonvec(n)=epsilon;

% Print
% fprintf([num2str(t*epsilon/E) ' ' num2str(epsilon) ' ' num2str(E) ' \n']);


end

% Plot Avg Energy Spectrum
Emean=Emean(31:78, :);
tvecint=tvec(31:78);
EmeanVec=zeros(500,1);
for n=1:500
	EmeanVec(n)=trapz(tvecint, Emean(:,n))/(tvecint(end)-tvecint(1));
end

	figure;
        title('Average Energy Spectrum');
        loglog(kline, EmeanVec);
        hold on;
        loglog(kline, 1.8191*114.85/34*1.4528*(kline.^4)./((1 + kline.^2).^(17/6)).*exp(-5.2.*kline.*0.0245));
        ylim([1e-18 1e2]);
	xlim([0.5, 300]);
	hold off;
        saveas(gcf,'meanE.fig');
        close;
        
% Plot Avg Enstrophy Spectrum
Omegamean=Omegamean(31:78, :);
OmegameanVec=zeros(500,1);
for n=1:500
	OmegameanVec(n)=trapz(tvecint, Omegamean(:,n))/(tvecint(end)-tvecint(1));
end

	figure;
        title('Average Enstrophy Spectrum');
        loglog(kline, OmegameanVec);
        hold on;
        loglog(kline, Omega0);
	ylim([1e-18 1e2]);
	xlim([0.5 300]);
        hold off;
        saveas(gcf,'meanOmega.fig');
        close;
% Plot TKE vs scaled time
figure;
title('TKE');
plot(tvec.*epsilonvec./Evec, Evec);
xlabel('t');
ylabel('TKE');
saveas(gcf,'TKE.fig');
close;

% Integral length
l=trapz(kline,(kline.^-1).*EmeanVec)/trapz(kline,EmeanVec);


save('mean.mat', 'EmeanVec', 'kline', 'OmegameanVec', 'tvec', 'Evec', 'epsilonvec', 'l');
