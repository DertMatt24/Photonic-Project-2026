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
        function draw_input(obj, in_ring)
            subplot(311);hold on;grid on;box on;plot(obj.time*1e9, in_ring,'k','LineWidth',2)
            xlabel('Time [ns]')
            ylabel('Input x(t)')
            xlim([obj.t_min obj.t_max])
            set(gca,'fontsize',12)
        end
        
        % This method plots the output signal over time.
        function draw_output(obj, out_ring, t, y)
            subplot(312);hold on;grid on;box on;
            plot(t*1e9, y,'k','LineWidth',2)
            grid on; plot(obj.time*1e9, out_ring,'r','LineWidth',2)
            xlabel('Time [ns]')
            ylabel('Output y(t)')
            xlim([obj.t_min obj.t_max])
            set(gca,'fontsize',12)
        end    
        
        % This method plots the output power over time
        function draw_power(obj, out_ring, t, y)
            subplot(313);hold on;grid on;box on;plot(t*1e9, abs(y).^2,'k','LineWidth',2)
            plot(obj.time*1e9, (out_ring).^2,'r','LineWidth',2)
            xlabel('Time [ns]')
            ylabel(' Output | y(t) |^2')
            xlim([obj.t_min obj.t_max])
            set(gca,'fontsize',12)
        end    
        
        % This method plots input, output signals and output power over
        % time.
        function input_output_power(obj, in_ring, out_ring, t, y)
            obj.draw_input(in_ring);
            obj.draw_output(out_ring, t, y);
            obj.draw_power(out_ring, t, y);      
        end

            
    end
    

    methods (Static)
        % This method plots the power spectrum in frequency domain.
        function spectrum_f(Df, H_drop, H_ODE, IN_ring)
            hold on; grid on, box on
            plot(Df/1e9,10*log10(abs(H_drop./max(abs(H_drop))).^2),'r','LineWidth',2)
            plot(Df/1e9,10*log10(abs(H_ODE./max(H_ODE)).^2),'b','LineWidth',2)
            plot(Df/1e9,10*log10(abs(IN_ring./max(abs(IN_ring))).^2),'k','LineWidth',2)
            set(gca,'fontsize',12)
            ylim([-30 0])
            xlim([-15 15])
            xlabel('Frequency [GHz]')
            ylabel('Spectrum [dB]')
        end
    end    
end    