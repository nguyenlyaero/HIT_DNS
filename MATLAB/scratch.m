[u,v,w] = generate_isotropic_turbulence(256,256,256,12345);

%%
[Eu,kline] = get_energy_spectrum(u,500);
[Ev,kline] = get_energy_spectrum(v,500);
[Ew,kline] = get_energy_spectrum(w,500);
figure;
loglog(kline, 1/2.*(Eu+Ev+Ew));
hold on;
%loglog(kline, 1.4528*(kline.^4)./((1 + kline.^2).^(17/6)).*exp(-5.2.*kline.*0.0245));

%%
kco=16*2/3;
uL=sharp_filter(u,kco);
uS=u-uL;
[ES,klineS] = get_energy_spectrum(uS,500);
[EL,klineL] = get_energy_spectrum(uL,500);

loglog(klineS, 3/2.*ES);
loglog(klineL, 3/2.*EL);

%%
dudx=ddx_hit(u); dudy=ddy_hit(u); dudz=ddz_hit(u);
dvdx=ddx_hit(v); dvdy=ddy_hit(v); dvdz=ddz_hit(v);
dwdx=ddx_hit(w); dwdy=ddy_hit(w); dwdz=ddz_hit(w);

epsilon=dudx.*dudx+dudy.*dudy+dudz.*dudz...
    +dvdx.*dvdx+dvdy.*dvdy+dvdz.*dvdz...
    +dwdx.*dwdx+dwdy.*dwdy+dwdz.*dwdz;
    

epsilon=1/114.85*mean(epsilon(:))