% Author: Vishnu Gunda / vsg0005@auburn.edu
% Date: 2024-10-01
% Assignment Name: hw02

classdef hw02
    methods (Static)

        function [c, n] = p1(f, a, b, epsilon, name, f_prime)
            % p1: Implement numerical methods to find the root of a function
            % :param f: function handle
            % :param a: real number, the left end of the interval
            % :param b: real number, the right end of the interval
            % :param epsilon: real number, the function tolerance
            % :param name: string, the name of the method
            % :param f_prime: function handle, the derivative of the function, only needed for Newton's method
            %
            % :return: 
            % c: real number, the root of the function
            % n: integer, the number of iterations

            if strcmp(name, 'bisection') % Bisection method
                n = 0;
                c = (a + b) / 2;
                while abs(f(c)) > epsilon
                    if f(a) * f(c) < 0
                        b = c;
                    else
                        a = c;
                    end
                    c = (a + b) / 2;
                    n = n + 1;
                end

            elseif strcmp(name, 'secant') % Secant method
                n = 0;
                c = b;  % initial guess c is the same as b
                x_prev = a;
                while abs(f(c)) > epsilon
                    x_new = c - f(c) * (c - x_prev) / (f(c) - f(x_prev));
                    x_prev = c;
                    c = x_new;
                    n = n + 1;
                end

            elseif strcmp(name, 'newton') % Newton-Raphson method
                n = 0;
                c = b;  % initial guess
                while abs(f(c)) > epsilon
                    fp = f_prime(c);  % Compute the derivative
                    if isnan(fp) || abs(fp) < 1e-12  % Check if the derivative is NaN or very small
                        fprintf('Newton-Raphson failed at iteration %d: f_prime(c) = %.5e\n', n, fp);
                        return;
                    end
                    c = c - f(c) / fp;
                    n = n + 1;
                end

            elseif strcmp(name, 'regula_falsi') % False position (Regula Falsi) method
                n = 0;
                c = a; % start with a
                while abs(f(c)) > epsilon
                    c = a - f(a) * (b - a) / (f(b) - f(a));
                    if f(a) * f(c) < 0
                        b = c;
                    else
                        a = c;
                    end
                    n = n + 1;
                end

            elseif strcmp(name, 'steffensen') % Steffensen's method
                n = 0;
                c = b; % initial guess
                while abs(f(c)) > epsilon
                    fc = f(c);  % Compute f(c)
                    if isnan(fc) || abs(fc) < 1e-12  % Check if f(c) is NaN or very small
                        fprintf('Steffensen failed at iteration %d: f(c) = %.5e\n', n, fc);
                        return;
                    end
                    fcf = f(c + fc);  % Compute f(c + f(c))
                    g_val = (fcf - fc) / fc;
                    if isnan(g_val) || abs(g_val) < 1e-12  % Check if g(c) is NaN or very small
                        fprintf('Steffensen failed at iteration %d: g(c) = %.5e\n', n, g_val);
                        return;
                    end
                    c = c - fc / g_val;
                    n = n + 1;
                end
            end
        end

        function p2()
            
        % summarize the iteration number for each method name in the table
        % Actual root is 1.37347
        
        %     |name          | iter | 
        %     |--------------|------|
        %     |bisection     |  32  |
        %     |secant        |  8   |
        %     |newton        |  6   |
        %     |regula_falsi  |  60  |
        %     |steffensen    |  145 |
            
        end

        function [c, n] = p3(f, a, b, epsilon)
            % Illinois method based on the provided algorithm in the image.
            % :param f: function handle
            % :param a: real number, the left end of the interval
            % :param b: real number, the right end of the interval
            % :param epsilon: real number, the function tolerance
            % :return: c: real number, the root of the function
            %          n: integer, the number of iterations


            % The Illinois Method Needed 9 Iterations

        
            % Initialization
            n = 0;          % Iteration counter
            max_iter = 1000; % Maximum number of iterations
            x0 = a;
            x1 = b;
            f0 = f(x0);     % f(a)
            f1 = f(x1);     % f(b)
            
            if f0 * f1 > 0
                error('The function must have opposite signs at a and b');
            end
            
            while true
                % Standard false position step
                x2 = x0 - f0 * (x1 - x0) / (f1 - f0);
                f2 = f(x2);
                n = n + 1;  % Increment iteration count
        
                % Check stopping criteria
                if abs(f2) < epsilon || abs(x2 - x1) < epsilon
                    c = x2;
                    return;
                end
                
                % If f1*f2 > 0, adjust until the sign changes (Illinois scaling)
                while f1 * f2 > 0
                    % Adjust f0 with λ = 1/2 (Illinois scaling)
                    f0 = f0 / 2;
                    x1 = x2;
                    f1 = f2;

                    % Update the variables according to the algorithm
                    x2=x0 - f0 * (x1 - x0)/(f1 - f0);
                    f2 = f(x2);
                    n = n + 1;
                    
                    % Check stopping criteria again
                    if abs(f2) < epsilon
                        c = x2;
                        return;
                    end
                end
                
                % If f1*f2 < 0, proceed with the false position step
                if f1 * f2 < 0
                    x0 = x1;
                    f0 = f1;
                    x1 = x2;
                    f1 = f2;
                end
                
                % Break out if maximum iterations exceeded
                if n >= max_iter
                    fprintf('Illinois method did not converge within %d iterations\n', max_iter);
                    c = NaN;
                    return;
                end
            end
        end

        function [c, n] = p4(f, a, b, epsilon)
            % Pegasus method for root finding
            % :param f: function handle
            % :param a: real number, the left end of the interval
            % :param b: real number, the right end of the interval
            % :param epsilon: real number, the function tolerance
            % :return: c: real number, the root of the function
            %          n: integer, the number of iterations


            % 52 Iterations were Needed for the Pegasus Method

        
            % Initialization
            n = 0;          % Iteration counter
            max_iter = 1000; % Maximum number of iterations
            x0 = a;
            x1 = b;
            f0 = f(x0);     % f(a)
            f1 = f(x1);     % f(b)
            
            % Ensure the function values at the boundaries have opposite signs
            if f0 * f1 > 0
                error('The function must have opposite signs at a and b');
            end
            
            while abs(x1 - x0) > epsilon && n < max_iter
                % Pegasus false position formula
                x2 = x0 - f0 * (x1 - x0) / (f1 - f0);
                f2 = f(x2);
                n = n + 1;
        
                % Stopping criterion
                if abs(f2) < epsilon
                    c = x2;
                    return;
                end
                
                % Update the interval
                if f0 * f2 < 0
                    % Root is in [x0, x2], so update x1
                    x1 = x2;
                    f1 = f2;
                else
                    % Root is in [x2, x1], apply Pegasus scaling
                    lambda = f1 / (f1 + f2);  % Pegasus method scaling
                    x0 = x2;
                    f0 = lambda * f0;  % Scale f0 using lambda
                end
                
                % Check for convergence
                if abs(x1 - x0) < epsilon
                    break;
                end
            end
            
            % Return the final root estimate
            c = (x0 + x1) / 2;
        
            % Check if the method reached the maximum iterations
            if n >= max_iter
                fprintf('Pegasus method did not converge within %d iterations\n', max_iter);
                c = NaN;
            end
        end
    end
end
