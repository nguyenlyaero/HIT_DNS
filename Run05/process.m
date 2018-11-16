clear; clc;
addpath('/scratch/06005/nguyenly/HIT_DNS/PadeOps_output');
addpath('/home1/06005/nguyenly/PadeOps/MATLAB');

% Read from Restart File
u0 = read_fortran_box('RESTART_Run05_u.000000', 512, 512, 512, 'double');
v0 = read_fortran_box('RESTART_Run05_v.000000', 512, 512, 512, 'double');
w0 = read_fortran_box('RESTART_Run05_w.000000', 512, 512, 512, 'double');

omega_x=ddy_hit(w0)-ddz_hit(v0);
omega_y=ddz_hit(u0)-ddx_hit(w0);
omega_z=ddx_hit(v0)-ddy_hit(u0);

[Omega_x,~]=get_energy_spectrum(omega_x,500);
[Omega_y,~]=get_energy_spectrum(omega_y,500);
[Omega_z,kline]=get_energy_spectrum(omega_z,500);

kline=kline.';

Omega0=(Omega_x+Omega_y+Omega_z).^(1/2);

n=66;
u = read_fortran_box(['Run05_uVel_t00' num2str(n,'%02d') '00.out'], 512, 512, 512, 'double');
v = read_fortran_box(['Run05_vVel_t00' num2str(n,'%02d') '00.out'], 512, 512, 512, 'double');
w = read_fortran_box(['Run05_wVel_t00' num2str(n,'%02d') '00.out'], 512, 512, 512, 'double');

tvec=t;

omega_x=ddy_hit(w)-ddz_hit(v);
omega_y=ddz_hit(u)-ddx_hit(w);
omega_z=ddx_hit(v)-ddy_hit(u);

[Omega_x,~]=get_energy_spectrum(omega_x,500);
[Omega_y,~]=get_energy_spectrum(omega_y,500);
[Omega_z,~]=get_energy_spectrum(omega_z,500);

Omegamean=(Omega_x+Omega_y+Omega_z).^(1/2);

[Eu,~] = get_energy_spectrum(u,500);
[Ev,~] = get_energy_spectrum(v,500);
[Ew,~] = get_energy_spectrum(w,500);

Emean=1/2*(Eu+Ev+Ew);

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

% Energy dissipation rate
dudx=ddx_hit(u); dudy=ddy_hit(u); dudz=ddz_hit(u);
dvdx=ddx_hit(v); dvdy=ddy_hit(v); dvdz=ddz_hit(v);
dwdx=ddx_hit(w); dwdy=ddy_hit(w); dwdz=ddz_hit(w);

epsilon=dudx.*dudx+dudy.*dudy+dudz.*dudz...
    +dvdx.*dvdx+dvdy.*dvdy+dvdz.*dvdz...
    +dwdx.*dwdx+dwdy.*dwdy+dwdz.*dwdz;
    
epsilon=1/114.85*mean(epsilon(:));


% Plot Avg Energy Spectrum

	figure;
        title('Average Energy Spectrum');
        loglog(kline, Emean);
        hold on;
        loglog(kline, 1.8191*114.85/34*1.4528*(kline.^4)./((1 + kline.^2).^(17/6)).*exp(-5.2.*kline.*0.0245));
        ylim([1e-18 1e2]);
	xlim([0.5 400]);
	hold off;
        saveas(gcf,'meanE.fig');
        close;
        
% Plot Avg Enstrophy Spectrum

	figure;
        title('Average Enstrophy Spectrum');
        loglog(kline, Omegamean);
        hold on;
        loglog(kline, Omega0);
	ylim([1e-18 1e2]);
	xlim([0.5 400]);
        hold off;
        saveas(gcf,'meanOmega.fig');
        close;

% Integral length
l=trapz(kline(2:end),(kline(2:end).^-1).*Emean(2:end))/trapz(kline(2:end),Emean(2:end));

EmeanVec=Emean; OmegameanVec=Omegamean;
save('mean.mat', 'EmeanVec', 'kline', 'OmegameanVec', 'l');
