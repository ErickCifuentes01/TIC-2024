clear all;
close all;
%% Create parameter sets for GFDM and OFDM
gfdm = get_defaultGFDM('TTI');
gfdm.K = 128;
gfdm.Kset = 0:33;  % Only allocate some subcarriers
gfdm.pulse = 'rrc_fd';
gfdm.a = 0.3;
gfdm.Mon =400;
gfdm.mu=10; %4QAM
nB = 1; % Number of GFDM blocks to generate
gfdm.L=0;%Number of overlapping subcarrier for frequency domain GFDM implementation..
gfdm.Kon=0;%Number of actives subcarriers.
gfdm.oQAM=1;
%% Generate the signals
% Currently only works without CP
%assert(~isfield(gfdm, 'Ncp') || gfdm.Ncp == 0);

% Allocate enough space for the signals
blockLen = gfdm.M*gfdm.K
fs=1.25e9;% frecuencia de muestreo
R=1;% impedancia de entrada


in_data=load("GFDM_RX\gfdm_signalrx_1024qam_18dbm.txt");
%in_data= load('GFDM.txt');
data=in_data';

%data = awgn(data,22);

% data=in_data;
% % %sincrononización señal-markers
ref_signal=zeros(length(data),1);  
xprbsr=PRBS([1 0 1 1 0 1 1 1 1],[9 8]);  % generador PRBS 7 para sincronismo
sync1=[0 xprbsr]; 
ref_signal(1:length(sync1)) = sync1';
rect_rx_data= data- mean(data);
Y1 = fft(rect_rx_data);
Y2 = fft(ref_signal);
Y = ifft(Y1.*conj(Y2));
[max1, nmax1] = max(abs(Y(1:(length(Y)/2))));
rx_data_sync = circshift(data, [-nmax1+1 0]);

 
data_gfdm=rx_data_sync(1:103296);
data_re=data_gfdm(1:51840);
data_im=data_gfdm(length(data_re)+1:length(data_gfdm));

rx_data1=data_re;
rx_data2=data_im;

%sincrononización real
xprbs=PRBS([1 0 1 1 0 1 1 1 1],[9 8]);  % generador PRBS 9 para sincronismo
sync1=[0 xprbs];             % generacion de la trama de sincronizacion[0+PBRS]
data_gfdm_real= rx_data1(length(sync1)+1:length(rx_data1));% 
 
% %sincrononización imaginaria
xprbsi=PRBS([0 1 0 1 0 1 1],[7 6]);  % generador PRBS 7 para sincronismo
sync12=[0 xprbsi];             % generacion de la trama de sincronizacion[0+PBRS]
data_gfdm_imag= rx_data2 (length(sync12)+1:length(rx_data2 ));% 
%  
% 
rx_data_gfdm=data_gfdm_real+data_gfdm_imag*1i;
xch=rx_data_gfdm;

in=do_demodulate(gfdm,xch);

 dhat = do_unmap(gfdm,in);

scatterplot( dhat)
 
 sh1 = do_qamdemodulate(dhat, gfdm.mu);
 
%simtx= load('simrandomicos.txt');
simtx = load("GFDM_TX\gfdm_bitstx_1024qam.txt");
 
  %% BER
 Bittx=de2bi( simtx,'left-msb');
 Bitrx=de2bi( sh1,'left-msb');
    bitx=Bittx(:);
    birx=Bitrx(:);
 errores= xor(bitx,birx);
 Nerror=sum(errores)
 Ntotal=length(bitx)
 BER=Nerror/(length(bitx)) 


EVM_RMS=comm.EVM;%Creación del objeto comm.EVM
EVM_RMS.ReferenceSignalSource="Estimated from reference constellation";%Definición del diagrama de constelación de referencia
EVM_RMS.ReferenceConstellation=qammod(0:2^gfdm.mu-1,2^gfdm.mu,'UnitAveragePower',true);%Diagrama de constelación de referencia con 2^BitsporSubportadora estados
EVM1=EVM_RMS(dhat);%Calculo del EVM de los símbolos de recepción
disp(['El EVM de la señal de recepcion es: ',num2str(EVM1)])

save("GFDM_RX\gfdm_bitsrx_1024qam_18dbm.txt","sh1","-ascii")
save GFDM_RX\gfdm_symrx_1024qam_18dbm.mat dhat





