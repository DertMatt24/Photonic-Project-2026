classdef mrr_asym
    %MRR_ASYM Class for simulating a Microring Resonator.
    
    properties
        R       % radius of the MRR
        neff    % effective index of the MRR waveguide
        ng      % group index of the MRR waveguide = nr1,nr2
        k_0     % zero power coupling coefficient at each single directional coupling point (Eq. 8, paper 4.6)
        A       % time scale factor
        alpha   % Loss factor

        Lb1     % Length of bus through port
        Lr1     % Length of ring arm through port
        Lb2     % Length of bus drop port
        Lr2     % Length of ring arm through port

        nb1     % waveguide group index of bus arm
        nb2     % waveguide group index of ring arm

        L_ring  % circumference of the ring [meters]
        tau     % round trip time 

        k1      % power coupling coefficient at through port
        k2      % power coupling coefficient at drop port
        tau_c   % cavity lifetime of the photons [seconds]
    end

    properties (Constant)
        c = 3e8; % light speed in vacuum [m/s]
    end
    
    methods
        % CLASS CONSTRUCTOR
        function obj = mrr_asym(R_input, neff_in, ng_in, k0_in, Lb1_in, Lr1_in, Lb2_in, Lr2_in, nb1_in, nb2_in, A, loss_factor)
            if nargin == 12 % robustness, if no enough parameters are used the object is not built.
                obj.R = R_input;
                obj.neff = neff_in; % effective index of the MRR waveguide
                obj.ng = ng_in; % group index of the MRR waveguide
                obj.nb1 = nb1_in; % group index of the MRR waveguide
                obj.nb2 = nb2_in; % group index of the MRR waveguide

                obj.Lb1 = Lb1_in;
                obj.Lr1 = Lr1_in;
                obj.Lb2 = Lb2_in;
                obj.Lr2 = Lr2_in;

                obj.k_0 = k0_in;

                obj.A = A;
                obj.alpha = loss_factor;

                obj.L_ring = 2*pi*R_input;
                obj.tau = obj.L_ring/(obj.c/obj.neff);

                obj.k1 = [];
                obj.k2 = [];
                obj.tau_c = [];
            else
               error("mrr build failed, be sure to have inserted 4 parameters!") 
            end    
            
        end

        %% Paper 4.6, Eq. 8: Effective coupling coefficients of the interferometric coupler
        %   KAPPA Calculates the effective coupling coefficients of the interferometric coupler.
        %   [k1, k2] = KAPPA(obj, lambda0) computes the effective power coupling 
        %   coefficients for the two arms of an interferometric coupler based on 
        %   wavelength, waveguide dimensions, phase shifts, and propagation losses.
        %
        %   This implementation corresponds to Paper 4.6, Equation 8.
        %
        %   INPUTS:
        %       obj     - Object containing the device parameters:
        %                   * k_0: Initial coupling coefficient
        %                   * nb1, nb2: Effective indices of the bus arms
        %                   * ng: Group index of the ring arms
        %                   * Lb1, Lb2: Lengths of the bus arms
        %                   * Lr1, Lr2: Lengths of the ring arms
        %                   * alpha: Propagation loss coefficient (dB/m)
        %       lambda0 - Free-space wavelength (same units as lengths)
        %
        %   OUTPUTS:
        %       k1      - Effective coupling coefficient for arm 1
        %       k2      - Effective coupling coefficient for arm 2
        %
        %   REMARKS:
        %   Loss is converted from dB to Neper (Np) for physical field attenuation.

        function [k1, k2] = kappa(obj,lambda0)
 
            phi_b1 = (2 * pi * obj.nb1 * obj.Lb1) / lambda0; % phase shift along the bus arm 1
            phi_b2 = (2 * pi * obj.nb2 * obj.Lb2) / lambda0; % phase shift along the bus arm 2
            phi_r1 = (2 * pi * obj.ng * obj.Lr1) / lambda0; % phase shift along the ring arm 1
            phi_r2 = (2 * pi * obj.ng * obj.Lr2) / lambda0; % phase shift along the ring arm 2
            
            alpha_Np = obj.alpha * log(10) / 20;

            T_b1 = exp(-2 * alpha_Np * obj.Lb1); % power transmission factor, bus arm 1
            T_r1 = exp(-2 * alpha_Np * obj.Lr1); % power transmission factor, ring arm 1
            T_b2 = exp(-2 * alpha_Np * obj.Lb2); % power transmission factor, bus arm 2 
            T_r2 = exp(-2 * alpha_Np * obj.Lr2); % power transmission factor, ring arm 2
 
            k1 = obj.k_0 * (1 - obj.k_0) * (T_b1 + T_r1 + 2 * sqrt(T_b1 * T_r1) * cos(phi_b1 - phi_r1));
            k2 = obj.k_0 * (1 - obj.k_0) * (T_b2 + T_r2 + 2 * sqrt(T_b2 * T_r2) * cos(phi_b2 - phi_r2));
        end

        %% Paper 4.6 tuning a0, b0
        
        function q = Q_Wu(obj, f0, ng, etha)
            w0 = 2 * pi * f0;

            q = - w0 * ng * obj.L_ring / (obj.c * log(1-etha));
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
        
        %% through port function
        function h_through = h_through_f(obj, Df, delta_f)
            alpha_Np = obj.alpha * log(10) / 20; 
            beta = 2 * pi * (Df - delta_f) / obj.c * obj.neff;
            gamma = (alpha_Np + 1i*beta) * obj.L_ring;
        
            t1 = sqrt(1 - obj.k1);
            t2 = sqrt(1 - obj.k2);
        
            h_through = (t1 - t2*exp(-gamma)) ./ (1 - t1*t2*exp(-gamma));
        end

        %% Through port — versione LINEARIZZATA (Eq. 3/4 del paper)
        % Da chiamare con a0, b0 gia calcolati (via parameters_LTI).
        function H_ODE_through = h_ode_through(obj, Df, a0, b0, delta_f)
            delta_Df = Df - delta_f;
            H_ODE_through = (b0 + 1i*2*pi*delta_Df) ./ (a0 + 1i*2*pi*delta_Df);
        end
        
        
        %%
        % Computing the Free Spectral Range
        % ng: group index
        function fsr = FSR(obj, ng)
            fsr = obj.c / (ng * obj.L_ring);
        end
        
        % Finding the resonant frequency
        % neff: effective index
        % M: integer greater than 0
        function f0 = resonant_frequency(obj, neff, M)
            f0 = obj.c / (neff * obj.L_ring) * M;
        end    
        
        % Banwidth at -3dB from the peak
        % fsr: Free Spectral Range
        % inserted inside paper 4.6
        % DELTAomega3dB = 2*a0 --> DELTAf3dB = a0/pi
        function b3db = B3dB(obj, a0)
            b3db = a0 / pi; % in Hz
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

