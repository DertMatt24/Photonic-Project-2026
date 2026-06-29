classdef Phase_detuning_analysis
    % Class containing just static methods to implement phase detuning
    % analysis
    
    methods (Static)
        % Phase detuning analysis for rings that use drop port
        function analysis_drop(in_ring, Df, MRR, ng)
            % Preparing input
            IN_ring=fftshift(fft(in_ring));
            fsr = MRR.FSR(ng);
            b3db = MRR.B3dB(fsr);
            
            % frequency axis
            N_sweep = 1000; 
            % multiplying by b3db to have a scale on x-axis that scales with the microring
            delta_f_array = linspace(-20 * b3db, 20 * b3db, N_sweep); 

            % output vectors
            phase_detuning = zeros(1, N_sweep);
            P_out  = zeros(1, N_sweep);
            
            % considering input power as constant
            P_in = mean(abs(in_ring).^2);
            
            % for loop to compute output power at different frequencies
            for i = 1:N_sweep
                delta_f = delta_f_array(i);
                phase_detuning(i) = 2 * pi * (delta_f / fsr);
                
                % ignoring normalized result -> it was used just to plot
                % graphs in other scripts
                [H_drop, ~] = MRR.h_drop_f(Df, delta_f);
                Out_ring = IN_ring .* H_drop;
                out_ring = real(ifft(fftshift(Out_ring)));
                
                % computing mean power to avoid complex calculations
                P_out(i) = mean(abs(out_ring).^2);
            end
            
            % Power lost
            P_lost = P_in - P_out;
            
            % Plot graph to see power usage depending from phase detuning
            figure(1)
            graph_drawer.power_graph(phase_detuning, P_out, P_in);
            
            % Plot grph to see power loss depending from phase detuning
            figure(2)
            graph_drawer.power_loss(phase_detuning, P_lost, P_in);
        end
        
        % Phase detuning analysis for MRR that use through port
        function analysis_through(in_ring, Df, MRR, ng)
            % Preparing input
            IN_ring=fftshift(fft(in_ring));
            fsr = MRR.FSR(ng);
            b3db = MRR.B3dB(fsr);
            
            % frequency axis
            N_sweep = 1000; 
            delta_f_array = linspace(-20 * b3db, 20 * b3db, N_sweep);

            % output vectors
            phase_detuning = zeros(1, N_sweep);
            P_out  = zeros(1, N_sweep);
            
            % considering input power as constant
            P_in = mean(abs(in_ring).^2);
            
            % for loop to compute output power at different frequencies
            for i = 1:N_sweep
                delta_f = delta_f_array(i);
                phase_detuning(i) = 2 * pi * (delta_f / fsr);
                
                % ignoring normalized result -> it was used just to plot
                % graphs in other scripts
                [H_through, ~] = MRR.h_through_f(Df, delta_f);
                Out_ring = IN_ring .* H_through;
                out_ring = real(ifft(fftshift(Out_ring)));
                
                % computing mean power to avoid complex calculations
                P_out(i) = mean(abs(out_ring).^2);
            end
            
            % Power lost
            P_lost = P_in - P_out;
            
            % Plot graph to see power usage depending from phase detuning
            figure(1)
            graph_drawer.power_graph(phase_detuning, P_out, P_in);
            
            % Plot grph to see power loss depending from phase detuning
            figure(2)
            graph_drawer.power_loss(phase_detuning, P_lost, P_in);
        end
    end
end    