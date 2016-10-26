function [A_out, b_out]=estimate_stable_lds(data, options)
% ESTIMATE_STABLE_LDS fits a stable linear dynamical system x_dot = A*x + b
%   to data
%
%   ESTIMATE_STABLE_LDS solves the optimization problem
%    min (error^2)  subject to:  error == A*x + b - x_dot
%     A                          A'+A <= -options.eps_constraints*I
%                                options.eps_constraints*I <= P
%   
%   [A_OUT, b_out] = ESTIMATE_STABLE_LDS(DATA, OPTIONS) returns the system
%   matrix and the bias of a linear dynamical system
%   
%   INPUT PARAMETERS:
%   -data    data = [x; x_dot] and size(data) = [d*2,n_data_points], 
%            where d is the dimenstion of the input/output.
%   -options options.solver -- specifies the YALMIP solver.
%            options.eps_constraints -- specifies the eps for the
%                                       constraints
%            options.bias (true|false) -- specifies if the LDS is 
%                                         x_dot = A*x + b (true) or just 
%                                         x_dot = A*x (false)
%
%   OUTPUT PARAMETERS:
%   - A_out  estimated system matrix
%   - b_out  estimated bias
%
%
%   # Authors: Sina Mirrazavi and Jose Medina
%   # EPFL, LASA laboratory
%   # Email: jrmout@gmail.com


d=size(data,1)/2;
options_solver=sdpsettings('solver',options.solver);

% Solver variables
A = sdpvar(d,d,'full');
b = sdpvar(d,1);
error = sdpvar(d,size(data,2));
objective_function=sum((sum(error.^2)));

if options.bias == false
    b = zeros(d,1);
end

% Define constraints
C=[error == (A*data(1:d,:) + repmat(b,1,size(data,2)))-data(d+1:2*d,:)];
% Lyapunov LMI setting P=I
C = C + [A'+A <= -options.eps_constraints*eye(d,d)] ;

% Solve the optimization
sol = optimize(C,objective_function,options_solver);
if sol.problem~=0
    warning(sol.info);
end

A_out = value(A);
b_out = value(b);
