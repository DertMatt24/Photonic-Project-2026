classdef graph_drawer
    % This class draws the graphs needed in the main script

    properties
        t_min % minimum time represented in the graph
        t_max % maximum time represented in the graph
        time % vector of time
    end
    
    methods
        
        function obj = graph_drawer(t_min, t_max, time)
            obj.t_min = t_min;
            obj.t_max = t_max;
            obj.time = time;
        end
        
        % This method plots the input signal over time.
        function draw_input(obj, in_ring, A)
            subplot(311);hold on;grid on;box on;plot(obj.time*A, in_ring,'k','LineWidth',2)
            xlabel('Time [ns]')
            ylabel('Input x(t)')
            xlim([obj.t_min obj.t_max])
            set(gca,'fontsize',12)
        end
        
        % This method plots the output signal over time.
        function draw_output(obj, out_ring, t, y, A)
            subplot(312);hold on;grid on;box on;
            plot(t*A, y,'k','LineWidth',2)
            grid on; plot(obj.time*A, out_ring,'r','LineWidth',2)
            xlabel('Time [ns]')
            ylabel('Output y(t)')
            xlim([obj.t_min obj.t_max])
            set(gca,'fontsize',12)
        end    
        
        % This method plots the output power over time
        function draw_power(obj, out_ring, t, y, A)
            subplot(313);hold on;grid on;box on;plot(t*A, abs(y).^2,'k','LineWidth',2)
            plot(obj.time*A, (out_ring).^2,'r','LineWidth',2)
            xlabel('Time [ns]')
            ylabel(' Output | y(t) |^2')
            xlim([obj.t_min obj.t_max])
            set(gca,'fontsize',12)
        end    
        
        % This method plots input, output signals and output power over
        % time.
        function input_output_power(obj, in_ring, out_ring, t, y, A)
            obj.draw_input(in_ring, A);
            obj.draw_output(out_ring, t, y, A);
            obj.draw_power(out_ring, t, y, A);      
        end

            
    end
    

    methods (Static)
        % This method plots the power spectrum in frequency domain.
        function spectrum_f(Df, H_drop, H_ODE, IN_ring, A)
            hold on; grid on, box on
            plot(Df/A,10*log10(abs(H_drop./max(abs(H_drop))).^2),'r','LineWidth',2)
            plot(Df/A,10*log10(abs(H_ODE./max(H_ODE)).^2),'b','LineWidth',2)
            plot(Df/A,10*log10(abs(IN_ring./max(abs(IN_ring))).^2),'k','LineWidth',2)

            legend('MRR ODE', 'Ideal ODE', 'Input Signal Spectrum', ...
           'Location', 'southwest', 'FontSize', 10);

            set(gca,'fontsize',12)
            ylim([-30 0])
            xlim([-15 15])
            xlabel('Frequency [GHz]')
            ylabel('Spectrum [dB]')
        end

        function power_graph(phase_detuning, P_optical_out, P_in)
            % Converting output in dB
            P_out_dB = 10 * log10((P_optical_out ./ P_in));
            
            plot(phase_detuning, P_out_dB, 'b', 'LineWidth', 1.5);
            grid on;
            xlabel('Phase Detuning \Delta\phi [rad]');
            ylabel('Power transmitted [dB]');
            title('Output Power vs Phase Detuning');
            xlim([min(phase_detuning) max(phase_detuning)]);
            ylim([max(P_out_dB)-50, max(P_out_dB)+2]); 
        end

        function power_loss(phase_detuning, P_lost, P_in)
            % Converting power lost in dB
            P_lost_dB = 10 * log10((P_lost ./ P_in));
       
            plot(phase_detuning, P_lost_dB, 'r', 'LineWidth', 2, 'DisplayName', 'Power loss [dB]');
            hold on; 
        
            % finding minimum power lost
            [min_val, ~] = min(P_lost_dB);
            yline(min_val, '--k', sprintf('  min = %.2f dB', min_val), ...
                  'LabelVerticalAlignment', 'bottom', 'Alpha', 0.5, ...
                  'DisplayName', 'Minimum power used');

            grid on;
            xlabel('Phase Detuning \Delta\phi [rad]');
            ylabel('Power lost [dB]');
            title('Power loss vs Phase detuning');
            legend('show', 'Location', 'southwest');
            xlim([min(phase_detuning) max(phase_detuning)]);
            ylim([min(P_lost_dB)-2, max(P_lost_dB)+2]);
            
            hold off; 
        end
    end    
end    