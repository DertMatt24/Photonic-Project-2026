close all
clear all 
clc

k_eq = 62.5; % ns^-1
k = 62.5;  % ns^-1
A = 1e9;  %time scaling parameter  [nano 10^9]
A_time = A * 1e3;
% Input signal x(t)
C=4;

% ODE SOLVER: k = Bs -> A_time = A * 1e3;
% FWHM = 19.07 * 1e-3;
% INTEGRATOR: k >> Bs -> A_time = A * 1e4;
% FWHM = 19.07 * 1e-4;
% SCALER: k << Bs -> A_time = A;
% FWHM = 19.07;
FWHM = 19.07 * 1e-3; 

%x = Model_utils.step_function(C); % step function (Heaviside) 
%x = Model_utils.sin_function(A);
x = Model_utils.super_gaussian_function(FWHM,A);
%x = Model_utils.gaussian_chirped(19.07, A, C);
%x = Model_utils.arbitrary_signal(A);

%x = Model_utils.rectangular(2*1e-5/A, 1);

% Definition of the ODE
odefun = Model_utils.first_order_ode(A, k_eq, x);

% Initial condition
y0 = 0;

N = 1e5;
% Time span
interval = 100/A_time;
tspan = linspace(-interval , interval, N);

% Solve using ode45
[t, y] = ode45(odefun, tspan, y0);


% Implementing the ODE solver with a microring resonator

R=30e-6;  %radiurn of the MRR
neff=2.4;   %effective index of the MRR waveguide
MRR = mrr(R, neff, k, A, 0);

fsr = MRR.FSR(MRR.neff);
b3db = MRR.B3dB(fsr);

n_fsr_range = 1; 
dt_max = 1 / (2 * n_fsr_range * fsr);

dt = (max(t) - min(t)) / (N - 1);
if dt > dt_max
    N = ceil((max(t) - min(t)) / dt_max) + 1;
    time = linspace(min(t), max(t), N);
    dt = time(2) - time(1);
else
    time = linspace(min(t), max(t), N);
end

in_ring = x(time);
IN_ring=fftshift(fft(in_ring));

Df=linspace(-1/(2*dt),1/(2*dt),N);


%% computing output
delta_f = linspace(-20 * b3db, 20 * b3db, N);

delta_f = 0;
[H_drop, H_drop_norm] = MRR.h_drop_f(Df, delta_f);
H_ODE = MRR.h_ode(Df, delta_f);

Out_ring=IN_ring.*H_drop;
Out_ODE=IN_ring.*H_ODE;

Out_ring_plot = IN_ring.*H_drop_norm;

out_ring=real(ifft(ifftshift(Out_ring)));
out_ode=real(ifft(ifftshift(Out_ODE)));

out_ring_plot = real(ifft(ifftshift(Out_ring_plot)));

%% Computing power loss for each architectures
[p, db] = MRR.power_loss(in_ring, out_ring, dt)

%% generate plots
t_min=-1;t_max=20;

utils = graph_drawer(t_min, t_max, time);

figure(1)
utils.input_output_power(in_ring, out_ring_plot,t ,y, A_time);

figure(2)
graph_drawer.spectrum_fsr(Df, H_drop, H_ODE, IN_ring, fsr, n_fsr_range);
figure(3)
graph_drawer.spectrum_f(Df, H_drop, H_ODE, IN_ring, A_time);

% comuting rmse
sum_sq_diff = 0;
N_total = length(time); 
for i = 1:N:N_total
    end_idx = min(i + N - 1, N_total);
    
    chunk_y = interp1(t, y, time(i:end_idx), 'linear', 0);
    chunk_out = out_ring_plot(i:end_idx);
    
    sum_sq_diff = sum_sq_diff + sum((chunk_y(:) - chunk_out(:)).^2);
end

rmse = sqrt(sum_sq_diff / N_total);