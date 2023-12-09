clc;
clear;
close all;

% Filter parameters
a = 0.999;                  % Parameter 'a' of the given filter
prec = 8;                   % Precision of fractional part
N = 1e4;                    % Number of samples
M = 50;                     % Number of simulations

% Transfer function of filter B/A
B = 1;
A = [1 -a];

% Quantization step
delta = 2^(-prec);

ype = zeros(1, M);
for i = 1:M
    % Generate a random error signal
    xe = rand(1,N)*delta - delta/2;
    % Filter signal according to transfer function
    ye = filter(B,A,xe);
    % Calculate output power
    ype(i) = sumsqr(ye)/length(ye);
end

% Expected power
qp = (delta^2/12)*(1/(1-abs(a)^2));

% Plot simulated power against expected value
fig = figure;
hold on;
stem(1:1:M, qp*ones(1,M));
stem(1:1:M, ype);
xlabel('Simulation Number');
ylabel('Output Quantization Power (W)');
title('Simulation of Output Quantization Power');
legend('Theory', 'Simulation');
grid on;
saveas(fig, [pwd,'/q2_fig1'], 'epsc');
