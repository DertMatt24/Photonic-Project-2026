close all
clear all 


k = 0.5;  % ns^-1
A=1e9;  %time scaling parameter  [nano 10^9]

% Input signal x(t)
C=4;
%x = Model_utils.step_function(C); % step function (Heaviside) 
x = Model_utils.super_gaussian_function(45.8, 1);

% Definition of the ODE
odefun = Model_utils.first_order_ode(A, k, x);

% Initial condition
y0 = 1;

% Time span
tspan = [-100e-9 100e-9];

% Solve using ode45
[t, y] = ode45(odefun, tspan, y0);


% Implementing the ODE solver with a microring resonator

R=5000e-6;  %radiurn of the MRR
neff=1.5;   %effective index of the MRR waveguide
MRR = mrr(R, neff, k, A);

%numerical version of the input step function x(t)
N=1e5;
time=linspace(min(t),max(t),N);
dt=time(2)-time(1);
in_ring = x(time);
IN_ring=fftshift(fft(in_ring));

Df=linspace(-1/(2*dt),1/(2*dt),N);


%% computing output

H_drop = MRR.h_drop_f(Df, k);
H_ODE = MRR.h_ode(Df, k);

Out_ring=IN_ring.*H_drop;
Out_ODE=IN_ring.*H_ODE;

out_ring=ifft(fftshift(Out_ring));
out_oude=ifft(fftshift(Out_ODE));

%% Tunable k from paper 4.7 using voltage
% Table to link faster Voltage to k tuned values
% V = [0.0  , 0.9   , 1.0   , 1.1   , 1.3]; [V]
% k = [38.0 , 46.0  , 54.0  , 63.0  , 82.0]; [ns-1]

MRR_Yang = mrr(R, neff, k, A);
MRR_Yang.tuning_voltage(1.3, A);

H_drop_Yang = MRR_Yang.h_drop_f(Df, k);
H_ODE_Yang = MRR_Yang.h_ode(Df, k);

Out_ring_Yang=IN_ring.*H_drop_Yang;
Out_ODE_Yang=IN_ring.*H_ODE_Yang;

out_ring_Yang=ifft(fftshift(Out_ring_Yang));
out_oude_Yang=ifft(fftshift(Out_ODE_Yang));


%% generate plots
t_min=-1;t_max=20;

utils = graph_drawer(t_min, t_max, time);

figure(1)
utils.input_output_power(in_ring, out_ring,t ,y);

figure(2)
graph_drawer.spectrum_f(Df, H_drop, H_ODE, IN_ring);


%% generate plots Yang 4.7

figure(3)
utils.input_output_power(in_ring, out_ring_Yang,t ,y);

figure(4)
graph_drawer.spectrum_f(Df, H_drop_Yang, H_ODE_Yang, IN_ring);
