clear; clc;
% addpath('~/HIT_DNS/MATLAB/');
addpath('/scratch/06005/nguyenly/HIT_DNS/PadeOps_output');
addpath('/home1/06005/nguyenly/PadeOps/MATLAB');

Run=3;
N=256;
Re=114.85;
Nsample=28;
tvec=zeros(Nsample,1);
I1PdfMat=zeros(Nsample,500);
I2PdfMat=zeros(Nsample,500);

for n=1:Nsample
% Read from file
u = read_fortran_box(['Run0' num2str(Run, '%d') '_uVel_t00' num2str(23+2*(n-1),'%02d') '00.out'], N, N, N, 'double');
v = read_fortran_box(['Run0' num2str(Run, '%d') '_vVel_t00' num2str(23+2*(n-1),'%02d') '00.out'], N, N, N, 'double');
w = read_fortran_box(['Run0' num2str(Run, '%d') '_wVel_t00' num2str(23+2*(n-1),'%02d') '00.out'], N, N, N, 'double');
fid =fopen(['Run0' num2str(Run, '%d') '_info_t00' num2str(23+2*(n-1),'%02d') '00.out']);
t=fscanf(fid,'%f'); t=t(1);
tvec(n)=t;

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

ddx_u=ddx_hit(u); ddy_u=ddy_hit(u); ddz_u=ddz_hit(u);
ddx_v=ddx_hit(v); ddy_v=ddy_hit(v); ddz_v=ddz_hit(v);
ddx_w=ddx_hit(w); ddy_w=ddy_hit(w); ddz_w=ddz_hit(w);


% I1
I1=-(uS.*uS.*sharp_filter(ddx_u,kco) ...
    + vS.*vS.*sharp_filter(ddy_v,kco) ...
    + wS.*wS.*sharp_filter(ddz_w,kco) ...
    + uS.*vS.*sharp_filter(ddy_u+ddx_v,kco) ...
    + uS.*wS.*sharp_filter(ddz_u+ddx_w,kco) ...
    + vS.*wS.*sharp_filter(ddz_v+ddy_w,kco));

if n==1
    [I1Pdf,Iline1, bin1] = get_transfer_pdf(I1,500);
else
    I1Pdf=get_transfer_pdf(I1,bin1);
end

I1PdfMat(n,:)=I1Pdf;

% I2
I2=-(sharp_filter(uL.*uL,kco) - uL.*uL).*(ddx_u-sharp_filter(ddx_u,kco)) ...
    + (sharp_filter(vL.*vL,kco) - vL.*vL).*(ddy_v-sharp_filter(ddy_v,kco)) ...
    + (sharp_filter(wL.*wL,kco) - wL.*wL).*(ddz_w-sharp_filter(ddz_w,kco)) ...
    + (sharp_filter(uL.*vL,kco) - uL.*vL).*(ddy_u+ddx_v-sharp_filter(ddy_u+ddx_v,kco)) ...
    + (sharp_filter(uL.*wL,kco) - uL.*wL).*(ddz_u+ddx_w-sharp_filter(ddz_u+ddx_w,kco)) ...
    + (sharp_filter(vL.*wL,kco) - vL.*wL).*(ddz_v+ddy_w-sharp_filter(ddz_v+ddy_w,kco));

if n==1
    [I2Pdf,Iline2, bin2] = get_transfer_pdf(I2,500);
else
    I2Pdf=get_transfer_pdf(I2,bin2);
end

I2PdfMat(n,:)=I2Pdf;

fprintf('%d \n', n);
end

% Time Average
I1Pdf=zeros(500,1);
I2Pdf=zeros(500,1);
for i=1:500
    I1Pdf(i)=trapz(tvec, I1PdfMat(:,i))./(tvec(end)-tvec(1));
    I2Pdf(i)=trapz(tvec, I2PdfMat(:,i))./(tvec(end)-tvec(1));
end

figure;
semilogy(Iline1,I1Pdf);
ylim([1e-5; inf]);
xlabel('I_1');
ylabel('PDF');
title('I1');
saveas(gcf, 'I1Pdf.fig');

figure;
semilogy(Iline2,I2Pdf);
ylim([1e-5; inf]);
xlabel('I_2');
ylabel('PDF');
title('I2');
saveas(gcf, 'I2Pdf.fig');

save('energytransfer.mat', 'I1Pdf', 'I2Pdf', 'Iline1', 'Iline2');
