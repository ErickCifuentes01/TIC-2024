%% SIGNAL UFM
% © Codigo UFMC GIETEC Universidad Politécnica Salesiana
% Derechos Reservados. Prohibido la reproducción total o parcial
% Propietario: GRUPO DE INVESTIGACIÓN GIETEC (Universidad Politécnica
% Salesiana)
% Elaborado/Modificado: Berenice Arguero T.

%close all
%clear all
clc

% %%% PARÁMETROS DE TRANSMISIÓN
% 
% %Número de puntos FFT
% numFFT = 2^19;
% %Tamaño de las sub-bandas (>1)
% subbandSize = 100;
% %Número de las sub-bandas (numSubbands*subbandSize<=numFFT)
% numSubbands = 100;
% %Separación entre sub-bandas
% subbandOffset = numFFT/2-subbandSize*numSubbands/2;
% %Longitud del filtro
% filterLen = 43;
% %Atenuación del lóbulo lateral (dB)
% slobeAtten = 50;
% %Bits por subportadora M-QAM (2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM)
% bitsPerSubCarrier = 4;
% %Creación del filtro Dolph-Chebyshev
% prototypeFilter = chebwin(filterLen, slobeAtten);
% 
% %Factor de incremento para la amplitud de la señal
% amplitud = 50;
% %Vector de relleno de ceros
% zerovector = zeros(1,10e2);


% %%% RECUPERACIÓN DE LA SEÑAL
% 
% %Cargar la señal desde el archivo txt
% data = load("UFMC_RX\ufmc_rx_16qam_03db.txt");
% %Obtener el número de muestras
% samples = length(data);
% %Muestrear la señal obtenida
% SenalRx = data(1:1:length(data));

% %%% PROCESAMIENTO DE LA SEÑAL
%Asignar la señal recibida a una variable
SenalRx = ufmc_rx;
% %Operación FFT
% in_datafft=fft(in_data);
% %Retirar la componente DC
% in_datafft(1)=0;
% %Realizar la IFFT
% in_data=real(ifft(in_datafft));%Operación IFFT
%Retiro del relleno de ceros
in_data = SenalRx(length(zerovector)+1:length(SenalRx)-length(zerovector));
% %Señal de referencia para la longitud de la señal
% ref_signal=zeros(length(in_data),1);
%Crear la señal PBRS para el sincronismo
sPRBS=2*PRBS([0 1 0 1 0 1 1 ],[7 6]);
%Señal de sincronización
sync1=sPRBS;
% %Crear vector de longitud de la señal PBRS
% ref_signal(1:length(sync1)) = sync1';
% %Rectificar la señal
% rect_data= in_data- mean(in_data);
% %FFT de la señal rectificada
% Y1 = fft(rect_data);
% %FFT de la señal de referencia
% Y2 = fft(ref_signal);
% %Operación conjugada entre la FFT de las señales rectificada y referencia
% Y = ifft(Y1.*conj(Y2));
% %Obtener los máximos del resultante
% [max1, nmax1] = max(abs(Y(1:(length(Y)/2))));
% %Sincronizar la señal
% signal_sync = circshift(in_data, [-nmax1+1 0]);





%Señal UFMC en la recepción
%SenalUFMC_Rx=signal_sync;
SenalUFMC_Rx = in_data;
%Retirar la señal PBRS y de la señal de divisibilidad
SenalUFMC1= SenalUFMC_Rx(length(sync1)+1:length(SenalUFMC_Rx)-21); 
%Obtener la parte real de la señal recibida
UFMC_real=SenalUFMC1(1:length(SenalUFMC1)/2)/amplitud;
%Obtener la parte imaginaria de la señal recibidad
UFMC_imag=SenalUFMC1(length(UFMC_real)+1:length(SenalUFMC1))/amplitud;
%Construir la señal UFMC
Rx_UFMC=UFMC_real+UFMC_imag*1i;





%Excluir ventanas y filtro adicionales
PaddedRx = [Rx_UFMC; zeros(2*numFFT-numel(Rx_UFMC),1)];
%PaddedRx = Rx_UFMC;
%Realizar la FFT
SimbolosRx2 = fftshift(fft(PaddedRx));
%Seleccionar las muestras pares
SimbolosRx = SimbolosRx2(1:2:end);
%Seleccionar las subportadoras de datos.
SimbolosSubportadoras = SimbolosRx(subbandOffset+(1:numSubbands*subbandSize));

%Desempaquetado OFDM
RxFrec = [prototypeFilter.*exp(1i*2*pi*0.5*(0:filterLen-1)'/numFFT);zeros(numFFT-filterLen,1)];
%Ecualización zero-forcing
PrototipoFiltroFrec = fftshift(fft(RxFrec));
PrototipoFilteroInv = 1./PrototipoFiltroFrec(numFFT/2-subbandSize/2+(1:subbandSize));

%Ordenar los símbolos
SimbolosRxMat = reshape(SimbolosSubportadoras,subbandSize,numSubbands);%Ordenamiento de los símbolos de recepción
%Ecualizador por sub-banda y eliminación de la distorsión por filtro
EqualizadosMat = bsxfun(@times,SimbolosRxMat,PrototipoFilteroInv);
%Obtener símbolos ecualizados
SimbolosEqualizadosRx = EqualizadosMat(:);

%Demodulador M-QAM
qamDemod = comm.RectangularQAMDemodulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitOutput', true, ...
    'NormalizationMethod', 'Average power');

%Demodular los símbolos QAM
bits_rx = qamDemod(SimbolosEqualizadosRx);

%Guardar los bits de recepción
%save bitsrx_16qam_03db.txt bits_rx -ascii
%save('datosRx.mat')

%Cargar los bits de transmisión
%bits_tx = load("UMFC_TX\bitstx_16qam.txt");
bits_tx = load("bitstx_4qam.txt");
%Calcular los errores
errors = xor(bits_tx,bits_rx);
%Calcular el BER
BER = sum(errors)/length(bits_tx)

%eyediagram(bits_rx,2)
constDiagram = comm.ConstellationDiagram;

constDiagram(SimbolosEqualizadosRx)

