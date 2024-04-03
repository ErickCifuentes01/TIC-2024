% Copyright (c) 2014 TU Dresden
% All rights reserved.
% See accompanying license.txt for details.
%% This example shows how a random GFDM signal can be
% created. Additionally, an OFDM signal with an equal length is
% created, and the Spectrum of both is compared.

%% Create parameter sets for GFDM and OFDM
clear all;
close all;
setPath
gfdm = get_defaultGFDM('TTI');
gfdm.K =2^7;
gfdm.Kset = 0:33;  % Only allocate some subcarriers
gfdm.pulse = 'rrc_fd';
gfdm.a = 0.3;
gfdm.Mon = 400;
gfdm.mu=2; % orden de modulación 
nB =1; % Number of GFDM blocks to generate
gfdm.L=0;%Number of overlapping subcarrier for frequency domain GFDM implementation..
gfdm.Kon=0;%Number of actives subcarriers.
gfdm.oQAM=1;
%% Generate the signals
% Currently only works without CP
%assert(~isfield(gfdm, 'Ncp') || gfdm.Ncp == 0);

% Allocate enough space for the signals
blockLen = gfdm.M*gfdm.K
% sGFDM = zeros(nB * blockLen, 1);
% sOFDM = zeros(size(sGFDM));

for b = 1:nB
    s = get_random_symbols(gfdm);
   
   %save gfdm_bitstx_1024qam.txt s -ascii
   save simrandomicos.txt s -ascii
 
    A=do_qammodulate(s, gfdm.mu);
    %A= qammod(s,16);
    scatterplot(A)
%      Areal=real(A); save sim4QAMrealtx.txt Areal -ascii
%      Aimag=imag(A); save sim4QAMimagtx.txt Aimag -ascii
%     save A.mat
    D = do_map(gfdm, A);
          x = do_modulate(gfdm,D,'F');
      %sGFDM((b-1)*blockLen+(1:blockLen)) = x;
      
end
x1= x';
close all

  
x2=x;% señal a graficar
fs=1.25e9;% frecuencia de muestreo
R=1;% impedancia de entrada
[Spdx,f] = psd_signal(x2,fs,R); %grafica la Power Spectrum Density (PSD)

save gfdm_tx.mat x2




% figure
% [psd,f] = periodogram(x2, rectwin(length(x2)), (2^7)*2,1000, 'centered'); 
% plot(f,10*log10(psd));


gfdm_real=real(x);
gfdm_imag=imag(x);


%-------- Genero trama para sincronismo---------------------
xprbsr=PRBS([1 0 1 1 0 1 1 1 1],[9 8]);  % generador PRBS 7 para sincronismo
xprbsr =xprbsr * max(abs(x1));   % mas amplitud para diferenciarlos
gfdm_r=[0 xprbsr gfdm_real']; 
xprbsi=PRBS([0 1 0 1 0 1 1],[7 6]);  % generador PRBS 7 para sincronismo
xprbsri =xprbsi * max(abs(x1));   % mas amplitud para diferenciarlos
gfdm_i=[0 xprbsi gfdm_imag']; 

GFDM=[gfdm_r gfdm_i];

%save gfdm_signaltx_1024qam.txt GFDM -ascii 
save GFDM.txt GFDM -ascii

%concatenación  real e imaginario
GFDM=[zeros(1,10e3) gfdm_r gfdm_i zeros(1,10e3)]';
figure(4)
plot(GFDM);

%save gfdm_signaltx_awg_1024qam.txt GFDM -ascii       %guarda los datos PRBS+FOFDM a tx
% 
% 
% 
