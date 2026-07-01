close all
clear all 
clc

k_eq = 0.5; % ns^-1
k = 0.5;  % ns^-1
A = 1e9;  %time scaling parameter  [nano 10^9]

% Input signal x(t)
C=4;
x = Model_utils.step_function(C); % step function (Heaviside) 
%x = Model_utils.sin_function(A);
%x = Model_utils.super_gaussian_function(19.07, A);
%x = Model_utils.gaussian_chirped(19.07, A, C);
%x = Model_utils.arbitrary_signal(A);

% Definition of the ODE
odefun = Model_utils.first_order_ode(A, k_eq, x);

% Initial condition
y0 = 0;

% Time span
tspan = [-100e-9 100e-9];

% Solve using ode45
[t, y] = ode45(odefun, tspan, y0);


% Implementing the ODE solver with a microring resonator

R=5000e-6;  %radiurn of the MRR
neff=1.5;   %effective index of the MRR waveguide
MRR = mrr(R, neff, k, A, 0);


N=1e5;
time=linspace(min(t),max(t),N);
dt=time(2)-time(1);
in_ring = x(time);
IN_ring=fftshift(fft(in_ring));

Df=linspace(-1/(2*dt),1/(2*dt),N);


%% computing output
fsr = MRR.FSR(1.5);
b3db = MRR.B3dB(fsr);
delta_f = linspace(-20 * b3db, 20 * b3db, N);

delta_f = delta_f(1);
[H_drop, H_drop_norm] = MRR.h_drop_f(Df, delta_f);
H_ODE = MRR.h_ode(Df, delta_f);

Out_ring=IN_ring.*H_drop;
Out_ODE=IN_ring.*H_ODE;

Out_ring_plot = IN_ring.*H_drop_norm;

out_ring=real(ifft(fftshift(Out_ring)));
out_ode=real(ifft(fftshift(Out_ODE)));

out_ring_plot = real(ifft(fftshift(Out_ring_plot)));

%% Tunable k from paper 4.7 using voltage
% Table to link faster Voltage to k tuned values
% V = [0.0  , 0.9   , 1.0   , 1.1   , 1.3]; [V]
% k = [38.0 , 46.0  , 54.0  , 63.0  , 82.0]; [ns-1]

MRR_Yang = mrr(R, neff, k, A, 0);
MRR_Yang.tuning_voltage(0.0, A);

[H_drop_Yang, H_drop_Yang_norm] = MRR_Yang.h_drop_f(Df, 0);
H_ODE_Yang = MRR_Yang.h_ode(Df, 0);

Out_ring_Yang = IN_ring .* H_drop_Yang;
Out_ODE_Yang=IN_ring.*H_ODE_Yang;

Out_ring_Yang_plot = IN_ring .* H_drop_Yang_norm;

out_ring_Yang = real(ifft(ifftshift(Out_ring_Yang)));
out_ode_Yang=ifft(fftshift(Out_ODE_Yang));

out_ring_Yang_plot = real(ifft(fftshift(Out_ring_Yang_plot)));

%% Tunable k from paper 4.6 using heaters
dx = Model_utils.derivative(x);



%% Computing power loss for each architectures
[p, db] = MRR.power_loss(in_ring, out_ring, dt)
[p_Yang, db_Yang]= MRR.power_loss(in_ring, out_ring_Yang, dt)
% attraverso la through port per Wu

%% generate plots
t_min=-1;t_max=20;

utils = graph_drawer(t_min, t_max, time);

figure(1)
utils.input_output_power(in_ring, out_ring_plot,t ,y, A);

figure(2)
graph_drawer.spectrum_f(Df, H_drop, H_ODE, IN_ring, A);


%% generate plots Yang 4.7

figure(3)
utils.input_output_power(in_ring, out_ring_Yang_plot,t ,y, A);

figure(4)
graph_drawer.spectrum_f(Df, H_drop_Yang, H_ODE_Yang, IN_ring, A);
