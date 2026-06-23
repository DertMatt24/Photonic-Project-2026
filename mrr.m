classdef mrr
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
    end
    
    properties (Constant)
        c = 3e8; % light speed in vacuum [m/s]
    end

    methods
        % CLASS CONSTRUCTOR
        % R_input: Radiurn of the MRR
        % neff_in: Effective index of the MRR waveguide
        % k: Inverse of the cavity life time (how much time a photon
        % remains trapped inside the ring)
        % A: Time scaling parameter used to normalize the ODE solver, for numerical stability.
        function obj = mrr(R_input, neff_in, k, A)
            if nargin == 4
                obj.R = R_input;
                obj.neff = neff_in; % effective index of the MRR waveguide

                obj.k_ring = k*A;
                obj.tau_c = 1/obj.k_ring; % cavity life time of the MRR

                obj.L_ring = 2*pi*R_input;

                obj.tau = obj.L_ring/(obj.c/obj.neff);
                obj.tau_n = obj.tau_c / obj.tau;

                obj.r = sqrt(obj.tau_n/(1+obj.tau_n));

            end    
            
        end
        
        % Computing the h_drop function in frequency domain
        function h_drop = h_drop_f(obj, Df, k)
            beta=2 * pi * Df / mrr.c * obj.neff;
            h_drop= 1/k*(1-obj.r^2)./(1-obj.r^2*exp(-1i*beta*obj.L_ring)); %frequency domain description of the MRR
        end    
        
        % Frequency domain descritpion of the ODE
        function H_ODE = h_ode(obj, Df, k)
            t_c = obj.tau_c;
            H_ODE=1/k*(1/t_c)./(1/t_c+1i*2*pi*Df);
        end    

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
        
        function quality_factor = Q(obj, f0)
            quality_factor = 2 * pi * f0 * obj.tau_c;
        end    

    end    

end
