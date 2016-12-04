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
%            options.attractor       -- specifies a priori the atractor
%            options.weights         -- weighting factor for each sample
%
%   OUTPUT PARAMETERS:
%   - A_out  estimated system matrix
%   - b_out  estimated bias
%
%
%   # Authors: Jose Medina and Sina Mirrazavi
%   # EPFL, LASA laboratory
%   # Email: jrmout@gmail.com

% Default values
if ~isfield(options, 'solver')
    options.solver = 'sedumi';
end 
if ~isfield(options, 'verbose')
    options.verbose = 0;
end 
if ~isfield(options, 'warning')
    options.warning = 0;
end 
if ~isfield(options, 'eps_constraints')
    options.eps_constraints = 1e-3;
end 

options_solver=sdpsettings('solver',options.solver, ...
                           'verbose', options.verbose, ...
                           'warning', options.warning);

% Solver variables
d=size(data,1)/2;
A = sdpvar(d,d,'full');
b = sdpvar(d,1);
error = sdpvar(d,size(data,2));

% Objective
objective_function=sum((sum(error.^2)));
if isfield(options, 'weights')
    objective_function=sum(options.weights.*(sum(error.^2)));
end

% Constraints
C=[error == (A*data(1:d,:) + repmat(b,1,size(data,2)))-data(d+1:2*d,:)];
% Lyapunov LMI setting P=I -> Pos def nonsymmetric matrix
C = C + [A'+A <= -options.eps_constraints*eye(d,d)] ;

% Do not estimate the bias, set it to the one specified a priori
if isfield(options, 'attractor')
    if (size(options.attractor,1) ~= d)
        error(['The specified attractor should have size ' d 'x1']);
    else
        C = C + [A*options.attractor-b == zeros(d,1)];
    end
end

% Solve
sol = optimize(C,objective_function,options_solver);
if sol.problem~=0
    warning(sol.info);
end

A_out = value(A);
b_out = value(b);
