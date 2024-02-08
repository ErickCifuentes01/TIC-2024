%Funcion para graficar la densidad de potencia espectral de una señal
%x muestras en tiempo
%fs frecuencia de muestreo en segundos
%R resistencia por defecto use 1 ohm para normalizar
function [Spdx,f] = psd_signal(x,fs,R)

f=-fs/2:fs/(length(x)-1):fs/2;
Z=(fftshift(fft(x))/(length(x)));
Spdx=10*log10((abs(Z).^2)/R*1000);
figure;
plot(f,Spdx);
title('Densidad de Potencia Espectral')
xlabel('Frecuencia (Hz)')
ylabel('Potencia (dBm)')