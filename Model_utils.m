classdef Model_utils
    
    %% INPUT SIGNAL GENERATOR
    methods (Static)
        % Step function at high C starting at time t0, t0 is 0 by default
        function x = step_function(C, t0)
            arguments
                C (1,1) double          
                t0 (1,1) double = 0  % optional, default 0
            end
            x = @(t) C * (t > t0);
        end
        
        % Sin function
        function x = sin_function(A)
            x = @(t) sin(A*t); %sinusoidal signal
        end    
           
        
        % Super gaussian function
        function x = super_gaussian_function(FWHM, A, m)
            arguments
                % FWHM in paper 4.7: FWHM = 19.07 ps m = 1 gaussian
                % FWHM in paper 4.7: FWHM = 41.54 ps for supergaussian
                FWHM (1,1) double
                A (1,1) double % time scaling parameter
                m (1,1) double = 1 % by dedfault is a gaussian
            end    
            FWHM = FWHM/A;
            x = @(t) exp( -log(2) .* ( (2 .* t) ./ FWHM ).^(2*m) );
        end    
        
        function x = gaussian_chirped(FWHM, A, C)
            FWHM = FWHM / A;
           
            sigma = FWHM / (2*sqrt(2*log(2)));
            
            x = @(t) exp(-(1 + 1j*C) .* t.^2 / (2*sigma^2));
        end    

        function x = arbitrary_signal(A)
            x = @(t) A.*t.*exp(-(A*t).^2).*cos(3*A*t);
        end    
        
        function xd = derivative(x)
            h = eps;
            xd = @(t) (x(t+h) - x(t-h))/ (2*h);
        end    
        
    end    
    
    %% ODE FUNCTION GENERATOR
    methods (Static)
        function odefun = first_order_ode(A, k, x)
            odefun = @(t,y) A*(x(t) - k*y);
        end

        function odefun = first_order_lti(a0, b0, x, x_d)
            odefun = @(t, y) -a0*y + x_d(t) + b0*x(t);
        end    
    end
    
    %% GENERAL UTILS
    methods (Static)
           
    end 
end    