clc;
clear;
close all;

% Filter parameters
h = [-0.375, 0.75, -0.375]; % Transfer function of filter
prec = 16;                  % Precision of fractional part
N = 1e4;                    % Number of samples

% Quantization step
delta = 2^(-prec);

% Generate a random error signal
xe = rand(1,N)*delta - delta/2;
% Filter signal according to transfer function
ye = filter(h,1,xe);

% Find frequency response of ye
w = -pi:2*pi/8192:pi;
[yp,wp] = freqz(ye,1,w);

% Smoothen the data for plotting
ypsd = smoothdata(10*log10(abs(yp)));

% Plot power spectral density in principal frequency range
fig = figure;
hold on;
plot(wp/pi, ypsd);
xlabel('Normalized Frequency (\times \pi rad/sec)');
ylabel('Power Spectral Density (dB)');
title('Output PSD of 16-Bit Quantization in a 3-Tap Filter');
grid on;
saveas(fig, [pwd,'/q3_fig1'], 'epsc');
