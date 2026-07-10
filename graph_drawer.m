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
            timescale = graph_drawer.find_time_label(A);
            subplot(311);hold on;grid on;box on;plot(obj.time*A, in_ring,'k','LineWidth',2)
            xlabel(timescale)
            ylabel('Input x(t)')
            xlim([obj.t_min obj.t_max])
            set(gca,'fontsize',12)
        end
        
        % This method plots the output signal over time.
        function draw_output(obj, out_ring, t, y, A)
            timescale = graph_drawer.find_time_label(A);

            subplot(312);hold on;grid on;box on;
            plot(t*A, y,'k--','LineWidth',2)
            grid on; plot(obj.time*A, out_ring,'r:','LineWidth',2)
            xlabel(timescale)
            ylabel('Output y(t)')
            xlim([obj.t_min obj.t_max])
            set(gca,'fontsize',12)
        end    
        
        % This method plots the output power over time
        function draw_power(obj, out_ring, t, y, A)
            timescale = graph_drawer.find_time_label(A);

            subplot(313);hold on;grid on;box on;plot(t*A, abs(y).^2,'k--','LineWidth',2)
            plot(obj.time*A, (out_ring).^2,'r:','LineWidth',2)
            xlabel(timescale)
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
            fscale = graph_drawer.find_frequency_label(A);
            
            hold on; grid on, box on
            plot(Df/A,10*log10(abs(H_drop./max(abs(H_drop))).^2),'r','LineWidth',2)
            plot(Df/A,10*log10(abs(H_ODE./max(H_ODE)).^2),'b','LineWidth',2)
            plot(Df/A,10*log10(abs(IN_ring./max(abs(IN_ring))).^2),'k','LineWidth',2)

            legend('MRR ODE', 'Ideal ODE', 'Input Signal Spectrum', ...
           'Location', 'southwest', 'FontSize', 10);

            set(gca,'fontsize',12)
            ylim([-30 0])
            xlim([-15 15])
            xlabel(fscale)
            ylabel('Spectrum [dB]')
        end
        
        % plotting the spectrum signal for fsr
        function spectrum_fsr(Df, H_drop, H_ODE, IN_ring, fsr, range)
            hold on; grid on, box on
            plot(Df/fsr,10*log10(abs(H_drop./max(abs(H_drop))).^2),'r','LineWidth',2)
            plot(Df/fsr,10*log10(abs(H_ODE./max(H_ODE)).^2),'b','LineWidth',2)
            plot(Df/fsr,10*log10(abs(IN_ring./max(abs(IN_ring))).^2),'k','LineWidth',2)

            legend('MRR ODE', 'Ideal ODE', 'Input Signal Spectrum', ...
           'Location', 'southwest', 'FontSize', 10);

            set(gca,'fontsize',12)
            ylim([-30 0])
            xlim([-range range])
            xlabel('Frequency/FSR')
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
            ylabel('Power loss [dB]');
            title('Power loss vs Phase detuning');
            legend('show', 'Location', 'southwest');
            xlim([min(phase_detuning) max(phase_detuning)]);
            ylim([min(P_lost_dB)-2, max(P_lost_dB)+2]);
            
            hold off; 
        end
        
        % Method to plot the power loss in function of k_values changing
        function plot_power_vs_k(k_values, power_k, k_perfect)
            arguments
                k_values (1,:) double 
                power_k  (1,:) double 
                k_perfect (1,1) double = NaN % value for which the MRR solves perfectly the input ODE
            end
            figure;
            plot(k_values, power_k, 'b-', 'LineWidth', 1.5);
            hold on;
            
            if ~ isnan(k_perfect)
                legend('Power loss', 'Location', 'best');
            else
                xline(min(mrr.k_Yang), '--k', 'k = 38 ns^-1', 'LabelVerticalAlignment', 'bottom', ...
                    'LineWidth', 1.5, 'DisplayName', 'Lower Limit (38 ns^-1)');
                
                xline(max(mrr.k_Yang), '--k', 'k = 82 ns^-1', 'LabelVerticalAlignment', 'bottom', ...
                    'LineWidth', 1.5, 'DisplayName', 'Upper Limit (82 ns^-1)');
                
                k_legend = 'k limits (Yang et al.)';
                legend('Power loss', k_legend, 'Location', 'best');
            end
            
            xlabel('k [ns^{-1}]');
            ylabel('Power loss [dB]');
            title('Power loss vs coefficient k');
            
            grid on;
            hold off;
        end
        
        % Method to plot rmse analysis
        function plot_rmse_analysis(k_values, rmse_values, k_perfect)
            arguments
                k_values (1,:) double 
                rmse_values  (1,:) double 
                k_perfect (1,1) double = NaN % value for which the MRR solves perfectly the input ODE
                                             % it corresponds to least rmse
                                             % value
            end

            figure;
            hold on;
            grid on;
            box on;
        
            plot(k_values, rmse_values, '-b', ...
                'LineWidth', 2, ...
                'MarkerEdgeColor', 'b', ...
                'MarkerFaceColor', 'w', ...
                'MarkerSize', 4);
        
            set(gca, 'FontSize', 12);
            xlim([min(k_values) max(k_values)]);
            
            xlabel('Coupling Coefficient k [ns^{-1}]', 'FontSize', 12);
            ylabel('RMSE', 'FontSize', 12);
            title('RMSE analysis on scaled output (1/k)', 'FontSize', 14, 'FontWeight', 'bold');
                
                
            if isnan(k_perfect)
                xline(min(mrr.k_Yang), '--k', 'k = 38 ns^-1', 'LabelVerticalAlignment', 'bottom', ...
                    'LineWidth', 1.5, 'DisplayName', 'Lower Limit (38 ns^-1)');
                
                xline(max(mrr.k_Yang), '--k', 'k = 82 ns^-1', 'LabelVerticalAlignment', 'bottom', ...
                    'LineWidth', 1.5, 'DisplayName', 'Upper Limit (82 ns^-1)');
                
                k_legend = 'k limits (Yang et al.)';
                legend('rmse', k_legend, 'Location', 'best');
            end
            
            

            hold off;
        end

        function timescale = find_time_label(A)
            scale = A / 1e9;
            
            if scale >= 1e6 
                timescale = 'Time [fs]';
            elseif scale >= 1e3
                timescale = 'Time [ps]';
            elseif scale >= 1
                timescale = 'Time [ns]';
            elseif scale >= 1e-3
                timescale = 'Time [us]';
            else
                timescale = 'Time [ms]';
            end
        end

        function fscale = find_frequency_label(A)
            scale = A / 1e9;
            
            if scale >= 1e6 
                fscale = 'Frequency [PHz]';
            elseif scale >= 1e3
                fscale = 'Frequency [THz]';
            elseif scale >= 1e2
                fscale = 'Frequency [100 GHz]';
            elseif scale > 10
                fscale = 'Frequency [10 GHz]';
            elseif scale >= 1
                fscale = 'Frequency [GHz]';
            elseif scale >= 1e-3
                fscale = 'Frequency [MHz]';
            else
                fscale = 'Frequency [KHz]';
            end
        end

    end    
end    