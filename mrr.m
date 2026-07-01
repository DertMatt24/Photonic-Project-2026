classdef mrr < handle
    % MRR Class for simulating a Microring Resonator.

    properties
        R       % radius of the MRR
        neff    % effective index of the MRR waveguide 
        k_ring  % scaled coupling coefficient (linked to the ODE model)
        tau_c   % cavity lifetime of the photons [seconds]
        L_ring  % circumference of the ring [meters]
        tau     % round trip time 
        tau_n   % dimensionless ratio between cavity lifetime and round trip time
        r       % coupling coefficient of the directional coupler of the MRR
        A       % time scale factor
        alpha   % Loss factor
        k_0     % power coupling coefficient at each single directional coupling point (Eq. 8, paper 4.6)
        
        % more on this on paper 4.7 Yang et al.
        k_v_fun % function that bonds voltage to k: k = f(V)
        inv_v_k % inverse of the aforementioned function: V = f^-1(k)
    end
    
    properties (Constant)
        c = 3e8; % light speed in vacuum [m/s]
        
        % data from paper 4.7, tuning k using these known values of k depending on V
        V_Yang = [0.0, 0.9, 1.0, 1.1, 1.3]; % [V]
        k_Yang = [38.0, 46.0, 54.0, 63.0, 82.0]; % [ns-1]
    end

    methods
        % CLASS CONSTRUCTOR
        % R_input: Radiurn of the MRR
        % neff_in: Effective index of the MRR waveguide
        % k: Inverse of the cavity life time (how much time a photon
        % remains trapped inside the ring)
        % A: Time scaling parameter used to normalize the ODE solver, for numerical stability.
        function obj = mrr(R_input, neff_in, k, A, loss_factor)
            if nargin == 5 % robustness, if no enough parameters are used the object is not built.
                obj.R = R_input;
                obj.neff = neff_in; % effective index of the MRR waveguide

                obj.k_ring = k*A;
                obj.tau_c = 1/obj.k_ring; % cavity life time of the MRR

                obj.L_ring = 2*pi*R_input;

                obj.tau = obj.L_ring/(obj.c/obj.neff);
                obj.tau_n = obj.tau_c / obj.tau;

                obj.r = sqrt(obj.tau_n/(1+obj.tau_n));

                obj.A = A;
                obj.alpha = loss_factor;

                % computing functions for paper 4.7
                [obj.k_v_fun, obj.inv_v_k] = mrr.k_v_function();
            else
               error("mrr build failed, be sure to have inserted 4 parameters!") 
            end    
            
        end
        
        % Function to tune the k parameter
        % k_new: new value tuned
        % A: time scale
        function obj = tuning_k(obj, k_new)
            % Changing k value and all other values that depends from k
            obj.k_ring = k_new * obj.A;
            obj.tau_c = 1/obj.k_ring;
            obj.tau_n = obj.tau_c / obj.tau;
            obj.r = sqrt(obj.tau_n/(1+obj.tau_n));
        end    
        
        %% Paper 4.7 tuning using voltage
        % Function to tune k using voltage as input
        % voltage: voltage value in [0.0, 0.9, 1.0, 1.1, 1.3]
        % A: time scale
        function obj = tuning_voltage(obj, voltage)
            if voltage < 0.0 || voltage > 1.3
                error("Voltage value not supported by current implementation. Voltage must be between 0.0 and 1.3 Volt")
            end

            k_new = obj.k_v_fun(voltage);
            obj.tuning_k(k_new);
                
        end    
         
        % This function simulates an external software that, given the
        % desired k we want to reach, it returns the voltage we need to
        % apply to the MRR.
        %
        % desired_k: k value we want to have on the MRR
        % 
        % v: voltage we have to apply to the MRR to reach the desired k
        function v = voltage_to_reach_k(obj, desired_k)
            v = obj.inv_v_k(desired_k);
        end

        %% Drop port
        % Computing the h_drop function in frequency domain
        function [h_drop, h_drop_normalized] = h_drop_f(obj, Df, delta_f)    
            k = obj.k_ring / obj.A;
            
            beta=2 * pi * (Df - delta_f) / mrr.c * obj.neff;
            gamma = (obj.alpha + 1i*beta) * obj.L_ring;

            h_drop= (1-obj.r^2)./(1-obj.r^2*exp(-gamma)); %frequency domain description of the MRR
            h_drop_normalized = 1/k * h_drop;
        end    
        
        % Frequency domain descritpion of the ODE
        function H_ODE_drop = h_ode(obj, Df, delta_f)
            k = obj.k_ring / obj.A; 
            t_c = obj.tau_c;
            H_ODE_drop=1/k*(1/t_c)./(1/t_c+1i*2*pi*(Df - delta_f));
        end    
        
        %% Through port
        % Through port transfer function
        function [h_through, h_through_norm] = h_through_f(obj, Df, delta_f)
            k = obj.k_ring / obj.A;
            beta = 2 * pi * (Df - delta_f) / mrr.c * obj.neff;
            gamma = (obj.alpha + 1i*beta) * obj.L_ring;

            h_through = obj.r * (1 - exp(-gamma)) ./ ...
                        (1 - obj.r^2 * exp(-gamma));
            h_through_norm = h_through * 1/k;
        end 
        
        function H_ODE_through = h_ode_through(obj, Df, a0, b0, delta_f)
            delta_Df = Df - delta_f;
            H_ODE_through = (b0 + 1i*2*pi*delta_Df) ./ (a0 + 1i*2*pi*delta_Df);
        end
        
        %%
        % Computing the Free Spectral Range
        % ng: group index
        function fsr = FSR(obj, ng)
            fsr = mrr.c / (ng*obj.L_ring);
        end
        
        % Finding the resonant frequency
        % neff: effective index
        % M: integer greater than 0
        function f0 = resonant_frequency(obj, neff, M)
            f0 = mrr.c / (neff * obj.L_ring) * M;
        end    
        
        % Banwidth at -3dB from the peak
        % fsr: Free Spectral Range
        function b3db = B3dB(obj, fsr)
            t_2 = 1 - obj.r^2;
            b3db = fsr * t_2 / pi;
        end

        % Computing the quality factor
        % f0: frequency [Hz]
        function quality_factor = Q(obj, f0)
            quality_factor = 2 * pi * f0 * obj.tau_c;
        end    
        
        %% Paper 4.6 tuning a0, b0
        
        function q = Q_Wu(obj, f0, ng, etha)
            w0 = 2 * pi * f0;

            q = - w0 * ng * obj.L_ring / (mrr.c * log(1-etha));
        end
        
        % This function computes the parameters definedd in the paper 4.6
        % Note: b1 is constant at one for first order LTI
        function [a0, b0] = parameters_LTI(obj, f0, Qi, Qe1, Qe2)
            w0 = 2*pi*f0;
            
            % doubling Q-factors to make the code readable
            d_Qi = Qi*2;
            d_Qe1 = Qe1*2;
            d_Qe2 = Qe2*2;
            
            % computing a0, b0 coefficients as the paper shows
            a0 = w0 * (1/d_Qi + 1/d_Qe1 + 1/d_Qe2);
            b0 = w0 * (1/d_Qi + 1/d_Qe2 - 1/d_Qe1);
        end

        %% Paper 4.6, Eq. 8: Effective coupling coefficients of the interferometric coupler
        %
        % Called once per coupler, calculate the related coupling
        % coefficient
        %
        % n_b : effective index of the bus arm
        % L_b : physical length of the bus arm
        % n_r : effective index of the ring arm
        % L_r :  physical length of the bus ring arm
        % lambda0  : wavelength at which the coefficient is evaluated
        %
        % Returns:
        % k        : effective power coupling coefficient of the coupler
        %
        function k = kappa(obj, n_b, L_b, n_r, L_r, lambda0)
            if isempty(obj.k_0)
                error("k_0 is not set. Define obj.k_0 (e.g. obj.kappa0 = 0.0441;) before calling kappa().")
            end
 
            phi_b = 2*pi*n_b*L_b/lambda0; % phase shift along the bus arm
            phi_r = 2*pi*n_r*L_r/lambda0; % phase shift along the ring arm
            
            alpha_Np = obj.alpha*log(10) / 20;

            T_b = exp(-2*alpha_Np*L_b); % power transmission factor, bus arm
            T_r = exp(-2*alpha_Np*L_r); % power transmission factor, ring arm
 
            k = obj.k_0*(1-obj.k_0) * (T_b + T_r + 2*sqrt(T_b*T_r)*cos(phi_b - phi_r));
        end


        %% Computing power difference
        
        % Computing power
        % x: signal received
        % dt: time domain where the signal is defined
        function power = power(obj, x, dt)
            power = sum(abs(x).^2)*dt;
        end    
        
        % Computing power loss in W and dB
        % x: input signal
        % y: output signal
        function [P_loss, P_loss_dB] = power_loss(obj, x, y, dt)
            P_in = obj.power(x, dt);
            P_out = obj.power(y, dt);

            P_loss = P_out - P_in;
            P_loss_dB = 10*log10(P_out/P_in);
        end
    end    
    
    methods (Static)
        function [fun, inv] = k_v_function()
            V_known = mrr.V_Yang;
            k_known = mrr.k_Yang;
            
            % from 0 to 0.9 we model the function as a second order polynom
            b_coeff = 5;     
            c_coeff = 38; % for V= 0.0 k is 38 ns-1
            a_coeff = (k_known(2) - c_coeff - b_coeff*V_known(2)) / V_known(2)^2;
            poly_part = @(v) a_coeff*v.^2 + b_coeff*v + c_coeff;
            
            % from 0.9 to 1.3 we connect the dots using linear interpolation
            V_linear = V_known(2:end);
            k_linear = k_known(2:end);
            
            % joining the two interpolation
            % [0.0, 0.9] -> polynomial 2nd order
            % (0.9, 1.3] -> linear
            fun = @(V) arrayfun(@(v) ...
                poly_part(v) * (v <= 0.9) + ...
                interp1(V_linear, k_linear, max(v, V_linear(1)), 'linear') * (v > 0.9), V);
            
            % calculating inverse function
            % pol inversion
            poly_inv = @(k) (b_coeff - sqrt(b_coeff^2 - 4*a_coeff*(c_coeff - k)))/ (-2*a_coeff); 
            % linear interpolation inversion
            linear_inv = @(k) interp1(k_linear, V_linear, max(k, k_linear(1)), 'linear');
            % computing inverse of all function
            inv = @(K) arrayfun(@(k) ...
                poly_inv(k) * (k <= 46) + ...
                linear_inv(k) * (k > 46), K);
        end
    end

        methods 
            function plot_k_v_function(obj)
                V_known = mrr.V_Yang;
                k_known = mrr.k_Yang;

                % Plot the function
                figure(1001); 
                hold on;
                V_dense = linspace(V_known(1), V_known(length(V_known)), 1000);

                plot(V_dense, obj.k_v_fun(V_dense), 'b-', 'LineWidth', 2, ...
                    'DisplayName', 'Function k(V)');
                scatter(V_known, k_known, 100, 'ro', 'filled', ...
                    'DisplayName', 'Known values');
                xlabel('Voltage [V]'); ylabel('k [ns^{-1}]');
                title('Function k(V)');
                legend('Location', 'northwest'); grid on;
                hold off;
                
                % Plot inverse
                figure(1002); 
                hold on;
                k_dense = linspace(k_known(1), k_known(length(k_known)), 1000);            
                
                plot(k_dense, obj.inv_v_k(k_dense), 'g-', 'LineWidth', 2, ...
                    'DisplayName', 'Inverse function V(k)');            
                scatter(k_known, V_known, 100, 'ro', 'filled', ...
                    'DisplayName', 'Known values');            
                xlabel('k [ns^{-1}]');
                ylabel('Voltage [Volt]');
                title('Inverse function');
                legend('Location', 'northwest');
                grid on;
                hold off;
            end
        end
end
