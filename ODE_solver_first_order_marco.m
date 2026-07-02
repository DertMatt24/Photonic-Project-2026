close all
clear all 
clc

A = 1e9;  %time scaling parameter  [nano 10^9]

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

lambda0 = 1550.391e-9; % resonance wavelenghts
f0 = mrr_asym.c / lambda0;

%% inizialization of MRR
MRR_Wu = mrr_asym(R_ring, neff, ng, k0, L_b1, L_r1, L_b2, L_r2, nb1, nb2, A, alpha);

% calcolo dei fattori k
[MRR_Wu.k1, MRR_Wu.k2] = MRR_Wu.kappa(lambda0);

Qi = MRR_Wu.Q_Wu(f0, ng, etha);       % Q dovuto alle perdite interne della cavità
Qe1 = MRR_Wu.Q_Wu(f0, ng, MRR_Wu.k1); % Q esterno dovuto alla Through Port
Qe2 = MRR_Wu.Q_Wu(f0, ng, MRR_Wu.k2); % Q esterno dovuto alla Drop Port

[a0, b0] = MRR_Wu.parameters_LTI(f0, Qi, Qe1, Qe2);

fsr = MRR_Wu.FSR(ng);
b3db = MRR_Wu.B3dB(a0);

% scaling the coefficients
a0_scaled = a0/A;
b0_scaled = b0/A;

% C=4;
%x = Model_utils.step_function(C); % step function (Heaviside) 
%x = Model_utils.sin_function(A);
%x = Model_utils.gaussian_chirped(19.07, A, C);
%x = Model_utils.arbitrary_signal(A);

% Input signal x(t)
x = Model_utils.super_gaussian_function(0.01, A);

% in our model we have a first derivative of the input
x_d = Model_utils.derivative(x);

% Definition of the ODE of paper 4.6
odefun = Model_utils.first_order_lti_scaled(a0, b0, x, x_d, A);

% Implementing the ODE solver with a microring resonator
% Initial condition
y0 = 0;

N = 1e5;
t_start = -100e-12;
t_end = 100e-12;

n_fsr_range = 4; 
dt_max = 1 / (2 * n_fsr_range * fsr);
dt_init = (t_end - t_start) / (N - 1);

if dt_init > dt_max
    N = ceil((t_end - t_start) / dt_max) + 1;
end

% Time span
tspan = linspace(t_start, t_end, N);
dt = tspan(2) - tspan(1);

tspan_ns = tspan * A;

% Solve using ode45 scaled to ns
[t_ns, y] = ode45(odefun, tspan_ns, y0);
t = t_ns / A;

in_ring = x(tspan);

IN_ring=fftshift(fft(in_ring));

Df = linspace(-1/(2*dt),1/(2*dt),N);

%% computing output
% delta_f = linspace(-20 * b3db, 20 * b3db, N);
% delta_f = delta_f(1);

delta_f = 0;

H_through_optical = MRR_Wu.h_through_f(Df,delta_f);
H_through_ODE = MRR_Wu.h_ode_through(Df, a0, b0, delta_f);

Out_through_optical_F = IN_ring .* H_through_optical;
Out_through_ODE_F = IN_ring .* H_through_ODE;

Out_through_optical = real(ifft(ifftshift(Out_through_optical_F)));
Out_through_ODE = real(ifft(ifftshift(Out_through_ODE_F)));

[p_Wu, db_Wu] = MRR_Wu.power_loss(in_ring, Out_through_optical, dt);
[p, db]= MRR_Wu.power_loss(in_ring, Out_through_ODE, dt);

fprintf('=== BILANCIO ENERGETICO ALLA THROUGH PORT ===\n');
fprintf('Variazione di Potenza (Modello Ottico): %.2f dB\n', db_Wu);
fprintf('Variazione di Potenza (Modello ODE - CMT):     %.2f dB\n', db);

%% generate plots
t_min=-0.05;t_max=0.05;

utils = graph_drawer(t_min, t_max, tspan);

figure(1)
utils.input_output_power(in_ring, Out_through_optical,t ,y, A);
% Perfezionamento dei titoli e legende per adattarsi al comportamento Through
subplot(311); title('Analisi Temporale alla Through Port (Modello Asimmetrico Wu)');
subplot(312); ylabel('Output y(t) & Through_{opt}'); legend('Soluzione ODE (y)', 'Modello Ottico Esatto');
subplot(313); ylabel('Output |y(t)|^2 & |Through_{opt}|^2'); legend('Potenza ODE |y|^2', 'Potenza Ottica Esatta');

figure(2)
graph_drawer.spectrum_f(Df, H_through_optical, H_through_ODE, IN_ring, A);
title('Spettro di Risposta in Frequenza alla Through Port');
legend('Filtro Ottico Esatto', 'Filtro Approssimato ODE', 'Spettro Segnale Ingresso');

% legenda AI slop perchè sono pigro
