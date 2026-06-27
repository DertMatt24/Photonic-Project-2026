classdef mrr < handle
    % MRR Class for simulating a Microring Resonator.

    properties
        R       % radiurn of the MRR
        neff    % effective index of the MRR waveguide 
        k_ring  % scaled coupling coefficient (linked to the ODE model)
        tau_c   % cavity lifetime of the photons [seconds]
        L_ring  % circumference of the ring [meters]
        tau     % round trip time 
        tau_n   % dimensionless ratio between cavity lifetime and round trip time
        r       % coupling coefficient of the directional coupler of the MRR
        A       % time scale factor
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
        function obj = mrr(R_input, neff_in, k, A)
            if nargin == 4 % robustness, if no enough parameters are used the object is not built.
                obj.R = R_input;
                obj.neff = neff_in; % effective index of the MRR waveguide

                obj.k_ring = k*A;
                obj.tau_c = 1/obj.k_ring; % cavity life time of the MRR

                obj.L_ring = 2*pi*R_input;

                obj.tau = obj.L_ring/(obj.c/obj.neff);
                obj.tau_n = obj.tau_c / obj.tau;

                obj.r = sqrt(obj.tau_n/(1+obj.tau_n));

                obj.A = A;
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
        
        % Function to tune k using voltage as input
        % voltage: voltage value in [0.0, 0.9, 1.0, 1.1, 1.3]
        % A: time scale
        function obj = tuning_voltage(obj, voltage, A)
            id = find(mrr.V_Yang == voltage);
            % if id is empty, it means that the voltage value inserted was
            % not documented inside the paper 4.7, so we return an error
            % due to unknown behaviour
            if ~isempty(id)
                obj.tuning_k(mrr.k_Yang(id)); % tuning to known k parameter
            else
                error("The voltage value is not supported. Supported values: [%s]", num2str(mrr.V))
            end    
                
        end    
        
        %% Drop port
        % Computing the h_drop function in frequency domain
        function [h_drop, h_drop_normalized] = h_drop_f(obj, Df)
            k = obj.k_ring / obj.A;
            beta=2 * pi * Df / mrr.c * obj.neff;
            h_drop= (1-obj.r^2)./(1-obj.r^2*exp(-1i*beta*obj.L_ring)); %frequency domain description of the MRR
            h_drop_normalized = 1/k * h_drop;
        end    
        
        % Frequency domain descritpion of the ODE
        function H_ODE_drop = h_ode(obj, Df)
            k = obj.k_ring / obj.A; 
            t_c = obj.tau_c;
            H_ODE_drop=1/k*(1/t_c)./(1/t_c+1i*2*pi*Df);
        end    
        
        %% Through port
        % Through port transfer function
        function [h_through, h_through_norm] = h_through_f(obj, Df)
            k = obj.k_ring / obj.A;
            beta = 2 * pi * Df / mrr.c * obj.neff;
            h_through = obj.r * (1 - exp(-1i*beta*obj.L_ring)) ./ ...
                        (1 - obj.r^2 * exp(-1i*beta*obj.L_ring));
            h_through_norm = h_through * 1/k;
        end 
        
        function H_ODE_through = h_ode_through(obj, Df, a0, b0)
            H_ODE_through = (b0 + 1i*2*pi*Df) ./ (a0 + 1i*2*pi*Df);
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
        
        % Computing the quality factor
        % f0: frequency [Hz]
        function quality_factor = Q(obj, f0)
            quality_factor = 2 * pi * f0 * obj.tau_c;
        end    
        
        %% Paper 4.6 tuning a0, b0
        
        function q = Q_Wu(obj, f0, ng, etha)
            w0 = 2 * pi * f0;

            q = w0 * ng * obj.L_ring / (mrr.c * log(1-etha));
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

end
