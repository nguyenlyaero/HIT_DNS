clear; clc;
% addpath('~/HIT_DNS/MATLAB/');
addpath('/scratch/06005/nguyenly/HIT_DNS/PadeOps_output');
addpath('/home1/06005/nguyenly/PadeOps/MATLAB');

Run=3;
N=256;
Re=114.85;
Nsample=28;
tvec=zeros(Nsample,1);
Term1Mat=zeros(Nsample, 500);
Term2Mat=zeros(Nsample, 500);
Term3Mat=zeros(Nsample, 500);
Term4Mat=zeros(Nsample, 500);
Term5Mat=zeros(Nsample, 500);
TermMat=zeros(Nsample, 500);
klineMat=zeros(Nsample, 500);

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

ddx_uL=ddx_hit(uL); ddy_uL=ddy_hit(uL); ddz_uL=ddz_hit(uL);
ddx_vL=ddx_hit(vL); ddy_vL=ddy_hit(vL); ddz_vL=ddz_hit(vL);
ddx_wL=ddx_hit(wL); ddy_wL=ddy_hit(wL); ddz_wL=ddz_hit(wL);

% Temporal term 1
Term1u=uL.*ddx_uS + vL.*ddy_uS + wL.*ddz_uS;
Term1u=Term1u - sharp_filter(Term1u, kco);

Term1v=uL.*ddx_vS + vL.*ddy_vS + wL.*ddz_vS;
Term1v=Term1v - sharp_filter(Term1v, kco);

Term1w=uL.*ddx_wS + vL.*ddy_wS + wL.*ddz_wS;
Term1w=Term1w - sharp_filter(Term1w, kco);

[Term1u, Term1v, Term1w]=project_divergence_free(Term1u, Term1v, Term1w);
[Term1uline, kline]=get_energy_spectrum(Term1u, 500);
[Term1vline, ~]=get_energy_spectrum(Term1v, 500);
[Term1wline, ~]=get_energy_spectrum(Term1w, 500);
Term1=Term1uline+Term1vline+Term1wline;
Term1Mat(n,:)=Term1;
klineMat(n,:)=kline;

% Temporal term 2
Term2u=uS.*ddx_uL + vS.*ddy_uL + wS.*ddz_uL;
Term2u=Term2u - sharp_filter(Term2u, kco);

Term2v=uS.*ddx_vL + vS.*ddy_vL + wS.*ddz_vL;
Term2v=Term2v - sharp_filter(Term2v, kco);


Term2w=uS.*ddx_wL + vS.*ddy_wL + wS.*ddz_wL;
Term2w=Term2w - sharp_filter(Term2w, kco);

[Term2u, Term2v, Term2w]=project_divergence_free(Term2u, Term2v, Term2w);
[Term2uline, ~]=get_energy_spectrum(Term2u, 500);
[Term2vline, ~]=get_energy_spectrum(Term2v, 500);
[Term2wline, ~]=get_energy_spectrum(Term2w, 500);
Term2=Term2uline+Term2vline+Term2wline;
Term2Mat(n,:)=Term2;

% Temporal term 3
Term3u=ddx_hit(uS.*uS - sharp_filter(uS.*uS, kco)) ...
    + ddy_hit(uS.*vS - sharp_filter(uS.*vS, kco)) ...
    + ddz_hit(uS.*wS - sharp_filter(uS.*wS, kco));

Term3v=ddx_hit(vS.*uS - sharp_filter(vS.*uS, kco)) ...
    + ddy_hit(vS.*vS - sharp_filter(vS.*vS, kco)) ...
    + ddz_hit(vS.*wS - sharp_filter(vS.*wS, kco));

Term3w=ddx_hit(wS.*uS - sharp_filter(wS.*uS, kco)) ...
    + ddy_hit(wS.*vS - sharp_filter(wS.*vS, kco)) ...
    + ddz_hit(wS.*wS - sharp_filter(wS.*wS, kco));

[Term3u, Term3v, Term3w]=project_divergence_free(Term3u, Term3v, Term3w);
[Term3uline, ~]=get_energy_spectrum(Term3u, 500);
[Term3vline, ~]=get_energy_spectrum(Term3v, 500);
[Term3wline, ~]=get_energy_spectrum(Term3w, 500);
Term3=Term3uline+Term3vline+Term3wline;
Term3Mat(n,:)=Term3;

% Temporal term 4
Term4u=ddx_hit(sharp_filter(uL.*uL, kco) - uL.*uL) ...
    + ddy_hit(sharp_filter(uL.*vL, kco) - uL.*vL) ...
    + ddz_hit(sharp_filter(uL.*wL, kco) - uL.*wL);

Term4v=ddx_hit(sharp_filter(vL.*uL, kco) - vL.*uL) ...
    + ddy_hit(sharp_filter(vL.*vL, kco) - vL.*vL) ...
    + ddz_hit(sharp_filter(vL.*wL, kco) - vL.*wL);

Term4w=ddx_hit(sharp_filter(wL.*uL, kco) - wL.*uL) ...
    + ddy_hit(sharp_filter(wL.*vL, kco) - wL.*vL) ...
    + ddz_hit(sharp_filter(wL.*wL, kco) - wL.*wL);

[Term4u, Term4v, Term4w]=project_divergence_free(Term4u, Term4v, Term4w);
[Term4uline, ~]=get_energy_spectrum(Term4u, 500);
[Term4vline, ~]=get_energy_spectrum(Term4v, 500);
[Term4wline, ~]=get_energy_spectrum(Term4w, 500);
Term4=Term4uline+Term4vline+Term4wline;
Term4Mat(n,:)=Term4;

% Temporal term 5
Term5u=(ddx_hit(ddx_uS) + ddy_hit(ddy_uS) + ddz_hit(ddz_uS))./Re;
Term5v=(ddx_hit(ddx_vS) + ddy_hit(ddy_vS) + ddz_hit(ddz_vS))./Re;
Term5w=(ddx_hit(ddx_wS) + ddy_hit(ddy_wS) + ddz_hit(ddz_wS))./Re;

[Term5uline, ~]=get_energy_spectrum(Term5u, 500);
[Term5vline, ~]=get_energy_spectrum(Term5v, 500);
[Term5wline, ~]=get_energy_spectrum(Term5w, 500);
Term5=Term5uline+Term5vline+Term5wline;
Term5Mat(n,:)=Term5;

% Temporal term
Termu=-Term1u - Term2u - Term3u + Term4u + Term5u;
Termv=-Term1v - Term2v - Term3v + Term4v + Term5v;
Termw=-Term1w - Term2w - Term3w + Term4w + Term5w;

[Termuline, ~]=get_energy_spectrum(Termu, 500);
[Termvline, ~]=get_energy_spectrum(Termv, 500);
[Termwline, ~]=get_energy_spectrum(Termw, 500);
Term=Termuline+Termvline+Termwline;
TermMat(n,:)=Term;

fprintf('%d \n', n);
end

kline=klineMat(1,:);

% Time Averaging
Term1=zeros(500,1);
Term2=zeros(500,1);
Term3=zeros(500,1);
Term4=zeros(500,1);
Term5=zeros(500,1);
Term=zeros(500,1);

for n=1:500
    Term1(n)=trapz(tvec, Term1Mat(:,n))./(tvec(end)-tvec(1));
    Term2(n)=trapz(tvec, Term2Mat(:,n))./(tvec(end)-tvec(1));
    Term3(n)=trapz(tvec, Term3Mat(:,n))./(tvec(end)-tvec(1));
    Term4(n)=trapz(tvec, Term4Mat(:,n))./(tvec(end)-tvec(1));
    Term5(n)=trapz(tvec, Term5Mat(:,n))./(tvec(end)-tvec(1));
    Term(n)=trapz(tvec, TermMat(:,n))./(tvec(end)-tvec(1));  
end

% Plot
figure;
loglog(kline, Term1);
hold on;
loglog(kline, Term2);
loglog(kline, Term3);
loglog(kline, Term4);
loglog(kline, Term, 'k--');
loglog(kline, Term5);
legend('Term 1 (sweeping)', 'Term 2 (straining)', 'Term 3 (nonlinear relaxation)', 'Term 4 (Leonard stress source/sink)', 'ddt (accelaration)', 'Term 5 (viscous stress)');
ylim([1e-7 1e2]);
hold off;
saveas(gcf, 'scalesplitting.fig');

save('scalesplitting.mat', 'Term1', 'Term2', 'Term3',  'Term4', 'Term5', 'Term', 'kline');
