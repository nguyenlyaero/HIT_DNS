clear; clc;
addpath('~/HIT_DNS/MATLAB/');
% addpath('/scratch/06005/nguyenly/HIT_DNS/PadeOps_output');
% addpath('/home1/06005/nguyenly/PadeOps/MATLAB');

Run=3;
N=256;
Re=114.85;


% Read from file
u = read_fortran_box(['Run0' num2str(Run, '%d') '_uVel_t00' num2str(78,'%02d') '00.out'], N, N, N, 'double');
v = read_fortran_box(['Run0' num2str(Run, '%d') '_vVel_t00' num2str(78,'%02d') '00.out'], N, N, N, 'double');
w = read_fortran_box(['Run0' num2str(Run, '%d') '_wVel_t00' num2str(78,'%02d') '00.out'], N, N, N, 'double');


% Scale splitting
kco=2/3*16;
uL=sharp_filter(u,kco);
uS=u-uL;
vL=sharp_filter(v,kco);
vS=v-vL;
wL=sharp_filter(w,kco);
wS=w-wL;

% Differentials
ddx_uS=ddx_hit(uS); ddy_uS=ddy_hit(uS); ddz_uS=ddz_hit(uS);
ddx_vS=ddx_hit(vS); ddy_vS=ddy_hit(vS); ddz_vS=ddz_hit(vS);
ddx_wS=ddx_hit(wS); ddy_wS=ddy_hit(wS); ddz_wS=ddz_hit(wS);

ddx_uL=ddx_hit(uL); ddy_uL=ddy_hit(uL); ddz_uL=ddz_hit(uL);
ddx_vL=ddx_hit(vL); ddy_vL=ddy_hit(vL); ddz_vL=ddz_hit(vL);
ddx_wL=ddx_hit(wL); ddy_wL=ddy_hit(wL); ddz_wL=ddz_hit(wL);

ddx_u=ddx_hit(u); ddy_u=ddy_hit(u); ddz_u=ddz_hit(u);
ddx_v=ddx_hit(v); ddy_v=ddy_hit(v); ddz_v=ddz_hit(v);
ddx_w=ddx_hit(w); ddy_w=ddy_hit(w); ddz_w=ddz_hit(w);


% I1
I1=-(uS.*uS.*ddx_uL ...
    + vS.*vS.*ddy_vL ...
    + wS.*wS.*ddz_wL ...
    + uS.*vS.*(ddy_uL+ddx_vL) ...
    + uS.*wS.*(ddz_uL+ddx_wL) ...
    + vS.*wS.*(ddz_vL+ddy_wL));

[I1Pdf,Iline1] = get_transfer_pdf(I1,500);

% I2
% I2=(sharp_filter(uL.*uL,kco) - uL.*uL).*ddx_uS ...
%     + (sharp_filter(vL.*vL,kco) - vL.*vL).*ddy_vS ...
%     + (sharp_filter(wL.*wL,kco) - wL.*wL).*ddz_wS ...
%     + (sharp_filter(uL.*vL,kco) - uL.*vL).*(ddy_uS+ddx_vS) ...
%     + (sharp_filter(uL.*wL,kco) - uL.*wL).*(ddz_uS+ddx_wS) ...
%     + (sharp_filter(vL.*wL,kco) - vL.*wL).*(ddz_vS+ddy_wS);

% Aditya
I2=(uL.*uL).*ddx_uS ...
    + (vL.*vL).*ddy_vS ...
    + (wL.*wL).*ddz_wS ...
    + (uL.*vL).*(ddy_uS+ddx_vS) ...
    + (uL.*wL).*(ddz_uS+ddx_wS) ...
    + (vL.*wL).*(ddz_vS+ddy_wS);

[I2Pdf,Iline2] = get_transfer_pdf(I2,500);


E=(u.*u+v.*v+w.*w)/2;
E=mean(E(:));

figure;
semilogy(Iline1,I1Pdf);
ylim([1e-5; inf]);
xlabel('I_1');
ylabel('PDF');
title('I1');

figure;
semilogy(Iline2,I2Pdf);
ylim([1e-5; inf]);
xlabel('I_2');
ylabel('PDF');
title('I2');


figure, surface(I1(:,:,45)','edgecolor','none'), daspect([1 1 1]), colorbar, title('I1');
figure, surface(I2(:,:,45)','edgecolor','none'), daspect([1 1 1]), colorbar, title('I2');

