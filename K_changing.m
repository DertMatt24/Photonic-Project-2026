close all
clear all 
clc

k_eq = 0.5; % ns^-1
k = 0.5;  % ns^-1
A=1e9;  %time scaling parameter  [nano 10^9]

% Input signal x(t)
C=4;
%x = Model_utils.step_function(C); % step function (Heaviside) 
x = Model_utils.sin_function(A);
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


%% creatinng MRR
R=5000e-6;  %radius of the MRR
neff=1.5;   %effective index of the MRR waveguide
MRR = mrr(R, neff, k, A, 0);

%% Computing input
N=1e5;
time=linspace(min(t),max(t),N);
dt=time(2)-time(1);
in_ring = x(time);
IN_ring=fftshift(fft(in_ring));

Df=linspace(-1/(2*dt),1/(2*dt),N);

%% ANALYSIS on CHANGING K_eq
% MRR extra parameters
k_perfect = k_eq; % k parameter of MRR to solve perfectly the ODE
is_tunable = false; % is the MRR tunable
port_used = "drop"; % the output port the MRR uses

            
N = 100;
%TODO: BALANCE THIS FUNCTION
k_values = linspace(max(k_perfect-50, eps), k_perfect + 50, N);

% preparing output vectors
% preparing power output vector
p_out_db = zeros(1, N);
% preparing rmse values
%rmse_values = zeros(1,N);
rmse_values_scaled = zeros(1,N);

% parameter useful for ODE
y0 = 0;
tspan = [-100e-9 100e-9];

% for loop changing k values
for i = 1 : length(k_values)

    if is_tunable
        if k_values(i) >= min(mrr.k_Yang) && k_values(i) <= max(mrr.k_Yang) 
            MRR.tuning_k(k_values(i));
        elseif k_values(i) <= min(mrr.k_Yang)
            MRR.tuning_k(min(mrr.k_Yang));
        elseif k_values(i) >= max(mrr.k_Yang)
            MRR.tuning_k(max(mrr.k_Yang));
        else
            error("Invalid k value reached")
        end
    end

    
    % suppose no phase detuning
    if strcmp(port_used, "drop")
        odefun = Model_utils.first_order_ode(A, k_values(i), x);
        [H_out, H_out_norm] = MRR.h_drop_f(Df, 0);
    elseif strcmp(port_used, "through")
        odefun = Model_utils.first_order_lti(a0,b0,x,x_d);
        [H_out, H_out_norm] = MRR.h_through_f(Df, 0);
    end
    
    sol = ode45(odefun, tspan, y0);
    

    Out_ring=IN_ring.*H_out;
    Out_ring_norm = IN_ring.*H_out_norm;
    
    out_ring=real(ifft(fftshift(Out_ring)));
    out_ring_norm = real(ifft(fftshift(Out_ring_norm)));
    
    [~, p_out_db(i)] = MRR.power_loss(in_ring, out_ring, dt);
    
    y = deval(sol, time); % making y length == out_ring_norm
    rmse_values_scaled(i) = sqrt(mean(abs(out_ring_norm - y).^2));
    %rmse_values(i) = sqrt(mean(abs(out_ring - y).^2));
end
% plotting results
if is_tunable 
    graph_drawer.plot_power_vs_k(k_values, p_out_db)
    before_commit.plot_rmse_analysis(k_values, rmse_values_scaled)
else
    graph_drawer.plot_power_vs_k(k_values, p_out_db, k_perfect)
    graph_drawer.plot_rmse_analysis(k_values, rmse_values_scaled, k_perfect)
end

