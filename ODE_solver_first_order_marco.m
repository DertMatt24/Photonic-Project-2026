close all
clear all 
clc

A = 1e9;  %time scaling parameter  [nano 10^9]

%% parameters definition from the paper:
L_ring = 178.98e-6; % circumference of the ring
R_ring = L_ring/(2*pi); % radius of the ring
ng = 4.1850;
neff = 2.45;

L_b1 = 116.8e-6;
L_r1 = 47.12e-6;
L_b2 = 47.12e-6;
L_r2 = 47.12e-6;

k0 = 0.0441;

%% heater providing power
eta1 = 6.35e-4;   % 1/mW, efficienza heater 1
eta2 = 6.47e-4;   % 1/mW, efficienza heater 2

% --> primo test: P_heater1 : 0, 1.52, 2.93, P_heater2 = 0.
% --> secondo test: P_heater1 : 0, P_heater2: 0, 2.92, 5,57
% --> terzo test: a0 const, b0 varied (0, 11.94), (1.60,10.28), (2.39,8.72)
% --> quarto test: a0 varied, b0 const

P_heater1_1 = 0; % mW <20
P_heater2_1 = 0; % mW <20

% secondo ordine
P_heater1_2 = 0; % mW <20
P_heater2_2 = 0; % mW <20

nb1_1 = ng + eta1 * P_heater1_1;
nb2_1 = ng + eta2 * P_heater2_1;

% secondo ordine
nb1_2 = ng + eta1 * P_heater1_2;
nb2_2 = ng + eta2 * P_heater2_2;

alpha = 800; % (db/m) loss factor
loss_dB_rt = alpha * L_ring; % loss in dB for round-trip
transmission_rt = 10^(-loss_dB_rt / 10); % survived power after round-trip 
etha = 1 - transmission_rt;

lambda0 = 1550.391e-9; % resonance wavelenghts
f0 = mrr_asym.c / lambda0;

%% inizialization of MRR
MRR_Wu = mrr_asym(R_ring, neff, ng, k0, L_b1, L_r1, L_b2, L_r2, nb1_1, nb2_1, A, alpha);

% calcolo dei fattori k
[MRR_Wu.k1, MRR_Wu.k2] = MRR_Wu.kappa(lambda0);

Qi = MRR_Wu.Q_Wu(f0, ng, etha);       % Q dovuto alle perdite interne della cavità
Qe1 = MRR_Wu.Q_Wu(f0, ng, MRR_Wu.k1); % Q esterno dovuto alla Through Port
Qe2 = MRR_Wu.Q_Wu(f0, ng, MRR_Wu.k2); % Q esterno dovuto alla Drop Port

[a0, b0] = MRR_Wu.parameters_LTI(f0, Qi, Qe1, Qe2)

MRR_Wu.tau_c = 1 / a0; % cavity lifetime 13 ps

fsr = MRR_Wu.FSR(ng);
b3db = MRR_Wu.B3dB(a0);

% scaling the coefficients
a0_scaled = a0/A;
b0_scaled = b0/A;

% C=4;
%x = Model_utils.step_function(2); % step function (Heaviside) 
%x = Model_utils.sin_function(A);
%x = Model_utils.gaussian_chirped(0.01907, A, 2);
%x = Model_utils.arbitrary_signal(A);

% Input signal x(t)
FWHM = 0.01907;
x = Model_utils.super_gaussian_function(FWHM, A);

% in our model we have a first derivative of the input
x_d = Model_utils.derivative(x, FWHM/A);

% Definition of the ODE of paper 4.6
odefun = Model_utils.first_order_lti_scaled(a0, b0, x, x_d, A);

% Implementing the ODE solver with a microring resonator
% Initial condition

y0 = 0;

N = 1e5;
t_start = -1e-9;
t_end = 1e-9;

fprintf('BS: %.3e >> tau_c: %.3e',t_end - t_start,MRR_Wu.tau_c);

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

y_norm = y ./ max(abs(y));

in_ring = x(tspan);

IN_ring=fftshift(fft(in_ring));

Df = linspace(-1/(2*dt),1/(2*dt),N);

%% computing output
% delta_f = linspace(-20 * b3db, 20 * b3db, N);
% delta_f = delta_f(1);

delta_f = 0;

H_through_optical = MRR_Wu.h_through_f(Df,delta_f);
H_through_ODE = MRR_Wu.h_ode_through(Df, a0, b0, delta_f);

% normalized

Out_through_optical_F = IN_ring .* H_through_optical;
Out_through_ODE_F = IN_ring .* H_through_ODE;

Out_through_optical = real(ifft(ifftshift(Out_through_optical_F)));
Out_through_ODE = real(ifft(ifftshift(Out_through_ODE_F)));

Out_through_optical_norm = Out_through_optical ./ max(abs(Out_through_optical));

[p_Wu, db_Wu] = MRR_Wu.power_loss(in_ring, Out_through_optical, dt);
[p, db]= MRR_Wu.power_loss(in_ring, Out_through_ODE, dt);

fprintf('=== BILANCIO ENERGETICO ALLA THROUGH PORT ===\n');
fprintf('Variazione di Potenza (Modello Ottico): %.2f dB\n', db_Wu);
fprintf('Variazione di Potenza (Modello ODE - CMT):     %.2f dB\n', db);

%% generate plots
t_min=-0.1;t_max=0.1;

utils = graph_drawer(t_min, t_max, tspan);

figure(1)
utils.input_output_power(in_ring, Out_through_optical_norm,t ,y_norm, A);
% Perfezionamento dei titoli e legende per adattarsi al comportamento Through
subplot(311); title('Temporal Analysis - Through Port (Asym model)');
subplot(312); ylabel('Output y(t) & Through_{opt}'); legend('Exact solution', 'Optical Power');
subplot(313); ylabel('Output |y(t)|^2 & |Through_{opt}|^2'); legend('Power ODE |y|^2', 'Optical Power');

% legenda AI slop perchè sono pigro

%% secondo ordine

MRR_Wu_2 = mrr_asym(R_ring, neff, ng, k0, L_b1, L_r1, L_b2, L_r2, nb1_2, nb2_2, A, alpha);

% calcolo dei fattori k
[MRR_Wu_2.k1, MRR_Wu_2.k2] = MRR_Wu_2.kappa(lambda0);

Qi_2 = MRR_Wu_2.Q_Wu(f0, ng, etha);       % Q dovuto alle perdite interne della cavità
Qe1_2 = MRR_Wu_2.Q_Wu(f0, ng, MRR_Wu_2.k1); % Q esterno dovuto alla Through Port
Qe2_2 = MRR_Wu_2.Q_Wu(f0, ng, MRR_Wu_2.k2); % Q esterno dovuto alla Drop Port

% riformulo i valori dati dal primo mrr
a10 = a0;
b10 = b0;

% calcolo i coefficienti del secondo mrr
[a20, b20] = MRR_Wu_2.parameters_LTI(f0, Qi_2, Qe1_2, Qe2_2)

% EDO del secondo ordine (sistema a due stati)
odefun2 = Model_utils.second_order_lti_scaled(a10, b10, a20, b20, x, x_d, A);
y0_2 = [0; 0];

[t_ns, y_2] = ode45(odefun2, tspan_ns, y0_2);

y1 = y_2(:,1); % uscita through primo mrr
y2 = y_2(:,2); % uscita through secondo mrr

% cross check frequency
H1_optical = MRR_Wu.h_through_f(Df, delta_f); 
H2_optical = MRR_Wu_2.h_through_f(Df, delta_f);
H_tot_optical = H1_optical .* H2_optical;

H1_ODE = MRR_Wu.h_ode_through(Df, a10, b10, delta_f);
H2_ODE = MRR_Wu_2.h_ode_through(Df, a20, b20, delta_f);
H_tot_ODE = H1_ODE .* H2_ODE;

Out2_optical = real(ifft(ifftshift(IN_ring .* H_tot_optical)));
Out2_ODE     = real(ifft(ifftshift(IN_ring .* H_tot_ODE)));

figure(2)
utils.input_output_power(in_ring, Out2_optical,t ,y2, A);
% Perfezionamento dei titoli e legende per adattarsi al comportamento Through
subplot(311); title('Temporal Analysis - Through Port (Asym model)');
subplot(312); ylabel('Output y(t) & Through_{opt}'); legend('ODE second order solution', 'Optical model');
subplot(313); ylabel('Output |y(t)|^2 & |Through_{opt}|^2'); legend('Power ODE second order |y|^2', 'Optical Power');


