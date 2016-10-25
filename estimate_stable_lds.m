function [A_out, b_out, P_out]=estimate_stable_lds(data, options)
% ESTIMATE_STABLE_LDS fits a stable linear dynamical system x_dot = A*x + b
%   to data
%
%   ESTIMATE_STABLE_LDS solves the optimization problem
%    min (error^2)  subject to:  error == A*x - x_dot
%     A                          A'*P+P*A <= -options.eps_constraints*I
%                                options.eps_constraints*I <= P
%   
%   [A_OUT, P_OUT] = ESTIMATE_STABLE_LDS(DATA, OPTIONS) offers two 
%   solutions specified in options.lyapunov_lmi:
%   - 'simplified': A simplified version where P = I. This option converges
%   faster and is suitable for simple and local dynamics
%
%   - 'full': the full optimization problem stated above which estimates
%   both A and P
%   
%   INPUT PARAMETERS:
%   -data    data = [x; x_dot] and size(data) = [d*2,n_data_points], 
%            where d is the dimenstion of the input/output.
%   -options options.solver -- specifies the YALMIP solver. 
%            options.lyapunov_lmi ('simplified'|'full')-- specifies 
%                                                        the lmi constraint
%            options.eps_constraints -- specifies the eps for the
%                                       constraints
%            options.bias (true|false) -- specifies if the LDS is 
%                                         x_dot = A*x + b (true) or just 
%                                         x_dot = A*x (false)
%
%   OUTPUT PARAMETERS:
%   - A_out  estimated system matrix
%   - P_out  estimated P matrix
%
%
%   # Authors: Sina Mirrazavi and Jose Medina
%   # EPFL, LASA laboratory
%   # Email: jrmout@gmail.com


d=size(data,1)/2;
options_solver=sdpsettings('solver',options.solver);

% Solver variables
A = sdpvar(d,d,'full');
P = sdpvar(d,d,'symmetric');
b = sdpvar(d,1);
error = sdpvar(d,size(data,2));
objective_function=sum((sum(error.^2)));

if options.bias == false
    b = zeros(d,1);
end

% Define constraints
C=[error == (A*data(1:d,:) + repmat(b,1,size(data,2)))-data(d+1:2*d,:)];
if strcmp(options.lyapunov_lmi,'simplified') == 1
    % P = I
    P = eye(d,d);
    % Lyapunov LMI
    C = C + [A'+A <= -options.eps_constraints*eye(d,d)] ;
elseif strcmp(options.lyapunov_lmi,'full') == 1 
    % Lyapunov LMI
    C = C + [A'*P+P*A <= -options.eps_constraints*eye(d,d)];
    % Positive definiteness of P
    C = C + [0.0001*eye(d,d)<=P];
else
    warning('Invalid option for the lyapunov lmi constraint.')
    return;
end

% Solve the optimization
sol = optimize(C,objective_function,options_solver);
if sol.problem~=0
    warn('An error occurred with the solver!')
end

A_out = value(A);
P_out = value(P);
b_out = value(b);
