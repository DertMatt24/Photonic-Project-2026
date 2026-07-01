close all
clear all 
clc

A = 1;  %time scaling parameter  [nano 10^9]

%% parameters definition from the paper:
L_ring = 178.98e-6; %circumference of the ring
R_ring = L_ring/(2*pi); %radius of the ring
ng = 4.1850;
neff = 2.45;

L_b1 = 116.8e-6;
L_r1 = 47.12e-6;
L_b2 = 47.12e-6;
L_r2 = 47.12e-6;

k0 = 0.0441;

% not explicitly defined inside the paper
nb1 = ng; % for now ng
nb2 = ng; % for now ng

alpha = 800; % (db/m) loss factor
loss_dB_rt = alpha * L_ring; % loss in dB for round-trip
transmission_rt = 10^(-loss_dB_rt / 10); % survived power after round-trip 
etha = 1 - transmission_rt;

lambda0 = 1550.391e-9;
f0 = mrr_asym.c / lambda0;

%% inizialization of MRR
MRR_Wu = mrr_asym(R_ring, neff, ng, k0, L_b1, L_r1, L_b2, L_r2, nb1, nb2, A, alpha);

% calcolo dei fattori k
[MRR_Wu.k1, MRR_Wu.k2] = MRR_Wu.kappa(lambda0);

Qi = MRR_Wu.Q_Wu(f0, ng, etha);       % Q dovuto alle perdite interne della cavità
Qe1 = MRR_Wu.Q_Wu(f0, ng, MRR_Wu.k1); % Q esterno dovuto alla Through Port
Qe2 = MRR_Wu.Q_Wu(f0, ng, MRR_Wu.k2); % Q esterno dovuto alla Drop Port

[a0, b0] = MRR_Wu.parameters_LTI(f0, Qi, Qe1, Qe2);
a0 = a0 / A; % scaled
b0 = b0 / A; % scaled

% Input signal x(t)
C=4;
x = Model_utils.step_function(C); % step function (Heaviside) 
%x = Model_utils.sin_function(A);
%x = Model_utils.super_gaussian_function(19.07, A);
%x = Model_utils.gaussian_chirped(19.07, A, C);
%x = Model_utils.arbitrary_signal(A);

% in our model we have a first derivative of the input
x_d = Model_utils.derivative(x);

% Definition of the ODE of paper 4.6
odefun = Model_utils.first_order_lti(a0, b0, x, x_d);

% Initial condition
y0 = 0;

% Time span
tspan = [-100e-9 100e-9];

% Solve using ode45
[t, y] = ode45(odefun, tspan, y0);



