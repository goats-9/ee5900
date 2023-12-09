clc;
clear;
close all;

% Filter parameters
a = 0.199;                  % Parameter 'a' of the given allpass filter
N = 256;                    % Number of samples computed
wl = [16];                  % Word length for fixed point implementations
prec = [2 4 8];             % Precision for each word length
w = -pi:(2*pi/8192):pi;     % Frequency range and step

% Array for legends
legend_arr = ["Floating point"];
for i = wl
    for j = prec
        str = append(num2str(i-j), "." ,num2str(j));
        legend_arr = [legend_arr, str];
    end
end

% For each implementation, we do the following:
% 1. Plot response with floating point precision
% 2. Plot response assuming 12,14,16-bit word length and
%    precision of 2,4,8 bits for each word length.

% Direct Form II
% --------------

% Floating Point Implementation

% Obtain time impulse response h[n]
y1 = ap_df2_float(a, N);

% Obtain frequency response H(jw)
[h1, w1] = freqz(y1, 1, w);

% Plot magnitude and phase response
fig = figure(1);
t = tiledlayout(2,1);
nexttile(1);
plot(w1/pi, 20*log10(abs(h1)));
hold on;
nexttile(2);
plot(w1/pi, angle(h1));
hold on;

% Fixed Point Implementation

for i = wl
    for j = prec
        % Obtain hq[n]
        y2 = ap_df2_fixed(a, i, j, N);
        % Find frequency response of quantized filter
        [h2, w2] = freqz(double(y2),1,w);
        % Plot magnitude and phase response
        nexttile(1);
        plot(w2/pi, 20*log10(abs(h2)));
        hold on;
        nexttile(2);
        plot(w2/pi, angle(h2));
        hold on;        
    end
end

% Further plotting commands

nexttile(1);
grid on;
xlabel('Normalized frequency (\times \pi rad/sec)');
ylabel('Magnitude (dB)');
legend(legend_arr);
title('Magnitude');
nexttile(2);
grid on;
xlabel('Normalized frequency (\times \pi rad/sec)');
ylabel('Magnitude (dB)');
legend(legend_arr);
title('Phase');

title(t, 'Direct Form II Allpass Digital Filter', 'fontweight', 'bold');
saveas(fig, [pwd,'/q1_fig1'], 'epsc');

% Reduced Multiplication
% ----------------------

% Floating Point Implementation

% Obtain time impulse response h[n]
y1 = ap_red_float(a, N);

% Obtain frequency response H(jw)
[h1, w1] = freqz(y1, 1, w);

% Plot magnitude and phase response
fig = figure(2);
t = tiledlayout(2,1);
nexttile(1);
plot(w1/pi, 20*log10(abs(h1)));
hold on;
nexttile(2);
plot(w1/pi, angle(h1));
hold on;

% Fixed Point Implementation

for i = wl
    for j = prec
        % Obtain hq[n]
        y2 = ap_red_fixed(a, i, j, N);
        % Find frequency response of quantized filter
        [h2, w2] = freqz(double(y2),1,w);
        % Plot magnitude and phase response
        nexttile(1);
        plot(w2/pi, 20*log10(abs(h2)));
        hold on;
        nexttile(2);
        plot(w2/pi, angle(h2));
        hold on;        
    end
end

% Further plotting commands

nexttile(1);
grid on;
xlabel('Normalized frequency (\times \pi rad/sec)');
ylabel('Magnitude (dB)');
legend(legend_arr);
title('Magnitude');
nexttile(2);
grid on;
xlabel('Normalized frequency (\times \pi rad/sec)');
ylabel('Magnitude (dB)');
legend(legend_arr);
title('Phase');

title(t, 'Reduced Allpass Digital Filter', 'fontweight', 'bold');
saveas(fig, [pwd,'/q1_fig2'], 'epsc');

% Floating point all pass filter, direct form II
function y = ap_df2_float(a, N)
    x = zeros(1,N);
    x(1) = 1;
    y = zeros(1, N);
    y(1) = -a;
    v = zeros(1, N);
    v(1) = 1;
    for n = 2:N
        v(n) = x(n) + a*v(n-1);
        y(n) = -a*v(n) + v(n-1);
    end
end

% Floating point all pass filter, reduced multiplications
function y = ap_red_float(a, N)
    x = zeros(1,N);
    x(1) = 1;
    y = zeros(1, N);
    y(1) = -a;
    for n = 2:N
        y(n) = a*(y(n-1)-x(n)) + x(n-1);
    end
end

% Fixed point all pass filter, direct form II
function y = ap_df2_fixed(a, wl, prec, N)
    x = zeros(1,N);
    x(1) = 1;
    y = fi(zeros(1, N), 1, wl, prec);
    y(1) = -a;
    v = fi(zeros(1, N), 1, wl, prec);
    v(1) = 1;
    for n = 2:N
        v(n) = x(n) + a*v(n-1);
        y(n) = -a*v(n) + v(n-1);
    end
end

% Fixed point all pass filter, reduced multiplications
function y = ap_red_fixed(a, wl, prec, N)
    x = zeros(1,N);
    x(1) = 1;
    y = fi(zeros(1, N), 1, wl, prec);
    y(1) = -a;
    for n = 2:N
        y(n) = a*(y(n-1)-x(n)) + x(n-1);
    end
end
