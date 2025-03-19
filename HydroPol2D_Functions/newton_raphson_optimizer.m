function [result, OF] = newton_raphson_optimizer(func, vars, initial_guess, max_iter, tol, damping)
% Newton-Raphson optimization for single or multiple variables in MATLAB.
%
% Parameters:
% func          - symbolic function to optimize
% vars          - symbolic variables as a vector, e.g., [x, y]
% initial_guess - initial point as a numeric vector
% max_iter      - maximum number of iterations (default: 100)
% tol           - tolerance for convergence (default: 1e-6)
% damping       - small factor to avoid singular Hessian issues (default: 1e-8)

if nargin < 4, max_iter = 100; end
if nargin < 5, tol = 1e-6; end
if nargin < 6, damping = 1e-8; end

% Calculate symbolic gradient and Hessian
grad_sym = gradient(func, vars);
hess_sym = hessian(func, vars);

% Convert symbolic functions to MATLAB functions
grad_func = matlabFunction(grad_sym, 'Vars', {vars});
hess_func = matlabFunction(hess_sym, 'Vars', {vars});
func_eval = matlabFunction(func, 'Vars', {vars});

x = initial_guess(:);  % Ensure column vector

for i = 1:max_iter

    % Evaluate gradient and Hessian
    grad_val = grad_func(x');
    grad_val = grad_val(:);  % Ensure column vector
    hess_val = hess_func(x');

    % Add damping
    if abs(det(hess_val)) < tol
        hess_val = hess_val + damping * eye(length(vars));
    end

    % Solve for update step
    dx = hess_val \ grad_val;

    % Update solution
    x_new = x - dx;

    % Check convergence
    if norm(x_new - x, 2) < tol
        fprintf('Converged in %d iterations.\n', i);
        result = x_new;
        OF = func_eval(x');        
        return;
    end

    x = x_new;
end

fprintf('Did not converge within the maximum number of iterations.\n');
result = x;
OF = func_eval(x');

end

% Example usage (uncomment to run):
% Define the symbolic objective function
% syms x_1 x_2;
% f = (1 - 8*x_1 + 7*x_1^2 - 7*x_1^3/3 + x_1^4/4) * x_2^2*exp(-x_2);
% 
% % Create the mesh of potential initial guesses for evaluating 
% xmin = [-10, -20];
% xmax = [30, 60];
% 
% xmed = (xmax - xmin)/2;
% n_tests = 200;
% 
% for i = 1:n_tests
%     x0_matrix(i,:) = rand(1,length(xmin)).*xmed;
% end
% 
% 
% for i = 1:n_tests   
%     [result(:,i),OF(i)] = newton_raphson_optimizer(f, [x_1, x_2], x0_matrix(i,:));
% end
% idx = find(OF == min(min(OF)),1,'first');
% disp('Optimal solution:');
% disp([result(:,idx); OF(idx)]);
