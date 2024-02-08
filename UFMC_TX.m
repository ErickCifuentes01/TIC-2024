%% SIGNAL UFM
% © Codigo UFMC GIETEC Universidad Politécnica Salesiana
% Derechos Reservados. Prohibido la reproducción total o parcial
% Propietario: GRUPO DE INVESTIGACIÓN GIETEC (Universidad Politécnica
% Salesiana)
% Elaborado/Modificado: Berenice Arguero T.

close all
clear all
clc

s = rng(211);       % Set RNG state for repeatability

%%% PARÁMETROS DE TRANSMISION

%Número de puntos FFT
numFFT = 2^14;
%Tamaño de las sub-bandas (>1)
subbandSize = 100;
%Número de las sub-bandas (numSubbands*subbandSize<=numFFT)
numSubbands = 100;
%Separación entre sub-bandas
subbandOffset = numFFT/2-subbandSize*numSubbands/2;
%Longitud del filtro
filterLen = 43;
%Atenuación del lóbulo lateral (dB)
slobeAtten = 50;
%Bits por subportadora M-QAM (2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM)
bitsPerSubCarrier = 2;
%Creación del filtro Dolph-Chebyshev
prototypeFilter = chebwin(filterLen, slobeAtten);
%Modulador M-QAM
qamMapper = comm.RectangularQAMModulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitInput', true, ...
    'NormalizationMethod', 'Average power');
%Factor de incremento para la amplitud de la señal
amplitud = 50;
%Creación del vector de datos
inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands);
%Creación de la señal UFMC
txSig = complex(zeros(numFFT+filterLen-1, 1));
%Vector para el registro de símbolos
sym_tx=[];
%Vector para los bits de transmisión
bitstx=[];

hFig = figure;
axis([-0.4 0.4 -80 0]);
hold on
grid on

xlabel('Frecuencia Normalizada');
ylabel('PSD (dBW/Hz)')
title(['UFMC, ' num2str(numSubbands) ' Subbandas x '  ...
    num2str(subbandSize) ' Subportadoras'])

%%% CONSTRUCCION DE LA SEÑAL UFMC

for bandIdx = 1:numSubbands
    %Bits de datos por sub-banda
    bitsIn = randi([0 1], bitsPerSubCarrier*subbandSize, 1);
    %Símbolos por sub-banda
    symbolsIn = qamMapper(bitsIn);
    %Concatenar bits en un solo vector
    inpData(:,bandIdx) = bitsIn;
    %Símbolos transmitidos
    sym_tx=[sym_tx; bitsIn];
    %Guardar los bits de transmisión
    save bitstx_4qam.txt sym_tx -ascii
    %Empaquetamiento de datos en un símbolo OFDM
    %Desplazamiento entre sub-bandas
    offset = subbandOffset+(bandIdx-1)*subbandSize;
    %Obtener los símbolos OFDM
    symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
        zeros(numFFT-offset-subbandSize, 1)];
    %Operación IFFT
    symifftOut = ifft(ifftshift(symbolsInOFDM));
    %Filtro para cada sub-banda desplazada en frecuencia
    bandFilter = prototypeFilter.*exp(1i*2*pi*(0:filterLen-1)'/numFFT* ...
        ((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2));
    %Señal con el filtro Dolph-Chevyshev
    filterOut = conv(bandFilter,symifftOut);
    % Plot power spectral density (PSD) per subband
    [psd,f] = periodogram(filterOut, rectwin(length(filterOut)), ...
                          numFFT*2, 1, 'centered'); 
    plot(f,10*log10(psd));
    %Señal UFMC que se va a transmitir
    txSig = txSig + filterOut;
end

hold off;

%Graficar diagrama de constelaciones de símbolos transmitidos
constDiagTx_UFMC = comm.ConstellationDiagram('ShowReferenceConstellation', ...
false, 'Position', figposition([20 60 20 25]), ...
'Title', 'UFMC Tx', ...
'Name', 'UFMC Tx', ...
'XLimits', [-1.5 1.5], 'YLimits', [-1.5 1.5]);
constDiagTx_UFMC(symbolsIn);

%%% SINCRONIZACION Y RELLENO DE CEROS

%Parte real de la señal UFMC
ufmc_real = real(txSig)';
%Parte imaginaria de la señal UFMC
ufmc_imag = imag(txSig)';
%Creación de la señal PBRS para el sincronismo
pbrs_signal = PRBS([0 1 0 1 0 1 1 ],[7 6]);
%Amplitud para diferenciar la señal PBRS
sync_signal = pbrs_signal*max(abs(txSig))*1.5;
%Garantizar la divisibilidad para 8
div8=zeros(21,1)';
%Concatenar la señal
ufmcout_tx = amplitud*[sync_signal ufmc_real ufmc_imag div8];
%Vector de relleno de ceros
zerovector = zeros(1,10e2);
%Señal UFMC con relleno de ceros
ufmc_tx = [zerovector ufmcout_tx zerovector]';
%Normalizar la señal UFMC
ufmc_tx = ufmc_tx/max(ufmc_tx);
%Guardar la señal UFMC
save ufmc_tx_4qam.txt ufmc_tx -ascii

%Crear señal de tiempo a 10Gbps
tbit = (1/(10e9))*1e6;
%Vector de tiempo en us
tiempo = linspace(0,tbit*length(ufmc_tx),length(ufmc_tx));

%Graficar la señal transmitida
figure
plot(tiempo,ufmc_tx);
title('Señal UFMC Transmitida - 4QAM');
xlabel('Tiempo (us)');
ylabel('Amplitud (V)');
grid on;





% 
% 
ufmc_rx = awgn(ufmc_tx,30);
ufmc_rx = ufmc_rx/max(ufmc_rx);

%Graficar la señal recibida
figure
plot(tiempo,ufmc_rx);
title('Señal UFMC Recibida - 4QAM');
xlabel('Tiempo (us)');
ylabel('Amplitud');
grid on;

%save ufmc_rx.txt ufmc_rx -ascii



% concatenación de señal
% paprf2=max((abs(ufmcout_tx)).^2)/mean((abs(ufmcout_tx)).^2);
% save signal.txt ufmcout_tx -ascii

%Graficar la señal en tiempo
% figure
% plot(ufmc_tx);
% scatterplot(txSig)
% 
% signalawgn = awgn(ufmcout_tx,10);
% 
% figure
% plot(signalawgn);
% 


rng(s);
