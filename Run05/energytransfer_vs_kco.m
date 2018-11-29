clear; clc;
% addpath('~/HIT_DNS/MATLAB/');
addpath('/scratch/06005/nguyenly/HIT_DNS/PadeOps_output');
addpath('/home1/06005/nguyenly/PadeOps/MATLAB');

Run=5;
N=512;
Re=114.85;
N0=48;
Nsample=10;
tvec=zeros(Nsample,1);
I1PdfMat=zeros(Nsample,500);
I2PdfMat=zeros(Nsample,500);
    
kcovec=2/3.*[8 16 32 64 128];
I1PdfMatvskco=zeros(length(kco), 500);
Iline1Matvskco=zeros(length(kco), 500);
I2PdfMatvskco=zeros(length(kco), 500);
Iline2Matvskco=zeros(length(kco), 500);
for j=1:length(kcovec)
    for n=1:Nsample
    % Read from file
    u = read_fortran_box(['Run0' num2str(Run, '%d') '_uVel_t00' num2str(N0+2*(n-1),'%02d') '00.out'], N, N, N, 'double');
    v = read_fortran_box(['Run0' num2str(Run, '%d') '_vVel_t00' num2str(N0+2*(n-1),'%02d') '00.out'], N, N, N, 'double');
    w = read_fortran_box(['Run0' num2str(Run, '%d') '_wVel_t00' num2str(N0+2*(n-1),'%02d') '00.out'], N, N, N, 'double');
    fid =fopen(['Run0' num2str(Run, '%d') '_info_t00' num2str(N0+2*(n-1),'%02d') '00.out']);
    t=fscanf(fid,'%f'); t=t(1);
    tvec(n)=t;

    % Scale splitting
    kco=kcovec(j);
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

    if n==1
        [I1Pdf,Iline1, bin1] = get_transfer_pdf(I1,500);
    else
        I1Pdf=get_transfer_pdf(I1,bin1);
    end

    I1PdfMat(n,:)=I1Pdf;

    % I2
    % I2=-(sharp_filter(uL.*uL,kco) - uL.*uL).*(ddx_u-sharp_filter(ddx_u,kco)) ...
    %     + (sharp_filter(vL.*vL,kco) - vL.*vL).*(ddy_v-sharp_filter(ddy_v,kco)) ...
    %     + (sharp_filter(wL.*wL,kco) - wL.*wL).*(ddz_w-sharp_filter(ddz_w,kco)) ...
    %     + (sharp_filter(uL.*vL,kco) - uL.*vL).*(ddy_u+ddx_v-sharp_filter(ddy_u+ddx_v,kco)) ...
    %     + (sharp_filter(uL.*wL,kco) - uL.*wL).*(ddz_u+ddx_w-sharp_filter(ddz_u+ddx_w,kco)) ...
    %     + (sharp_filter(vL.*wL,kco) - vL.*wL).*(ddz_v+ddy_w-sharp_filter(ddz_v+ddy_w,kco));

    % Aditya
    I2=(uL.*uL).*ddx_uS ...
        + (vL.*vL).*ddy_vS ...
        + (wL.*wL).*ddz_wS ...
        + (uL.*vL).*(ddy_uS+ddx_vS) ...
        + (uL.*wL).*(ddz_uS+ddx_wS) ...
        + (vL.*wL).*(ddz_vS+ddy_wS);

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
    title(['I1 (k_{co}=2/3*' num2str(2^(j+2),'%d') ')']);
    saveas(gcf, ['I1Pdf' num2str(2^(j+2),'%d') '.fig']);

    figure;
    semilogy(Iline2,I2Pdf);
    ylim([1e-5; inf]);
    xlabel('I_2');
    ylabel('PDF');
    title(['I2 (k_{co}=2/3*' num2str(2^(j+2),'%d') ')']);
    saveas(gcf, ['I2Pdf' num2str(2^(j+2),'%d') '.fig']);
    
    I1PdfMatvskco(j,:)=I1Pdf;
    Iline1Matvskco(j,:)=Iline1;
    I2PdfMatvskco(j,:)=I2Pdf;
    Iline2Matvskco(j,:)=Iline2;
end



save('energytransferVsKco.mat', 'I1PdfMatvskco', 'I2PdfMatvskco', 'Iline1Matvskco', 'Iline2Matvskco');
