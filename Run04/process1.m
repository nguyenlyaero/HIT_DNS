clear; clc;
addpath('/scratch/06005/nguyenly/HIT_DNS/PadeOps_output');
addpath('/home1/06005/nguyenly/PadeOps/MATLAB');

uMat=zeros(384, 384, 384, 37);
vMat=zeros(384, 384, 384, 37);
wMat=zeros(384, 384, 384, 37);
tVec=zeros(37,1);
for i=1:37
n=i+30;
u = read_fortran_box(['Run04_uVel_t00' num2str(n,'%02d') '00.out'], 384, 384, 384, 'double');
v = read_fortran_box(['Run04_vVel_t00' num2str(n,'%02d') '00.out'], 384, 384, 384, 'double');
w = read_fortran_box(['Run04_wVel_t00' num2str(n,'%02d') '00.out'], 384, 384, 384, 'double');
fid =fopen(['Run04_info_t00' num2str(n,'%02d') '00.out']);
t=fscanf(fid,'%f'); tVec(i)=t(1);

uMat(:,:,:,i)=u;
vMat(:,:,:,i)=v;
wMat(:,:,:,i)=w;
end

umean=zeros(384,384,384);
vmean=zeros(384,384,384);
wmean=zeros(384,384,384);
parfor i=1:384
    for j=1:384
        for k=1:384
            umean(i,j,k)=trapz(tVec,uMat(i,j,k,:))./(tVec(end)-tVec(1));
            vmean(i,j,k)=trapz(tVec,vMat(i,j,k,:))./(tVec(end)-tVec(1));
            wmean(i,j,k)=trapz(tVec,wMat(i,j,k,:))./(tVec(end)-tVec(1));
        end
    end
    
end

save('meanvel.mat', 'umean', 'vmean', 'wmean');





