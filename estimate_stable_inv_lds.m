function [A_inv, x_attractor]=estimate_stable_inv_lds(data, options, varargin)
% ESTIMATE_STABLE_LDS_INVERSE fits a stable linear dynamical system
% inverting the estimation problem. Instead of considering the standard
% model
%   x_dot = A*x + b = -A (x - x_attractor)
% it assumes that A is positive definite and performs the regression task
% in the inverse model, i.e.
%   x = x_attractor - A_inv x_dot
% 
% this model has the benefit of having the attractor as bias. As a drawback
% it performs worse than the standard regression procedure when it has very
% few data points. 
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
%            options.attractor         -- specifies a priori the atractor
%            options.weights           -- weighting factor for each sample
%
%   OUTPUT PARAMETERS:
%   - A_out  estimated system matrix
%   - b_out  estimated bias
%
%
%   # Author: Jose Medina
%   # EPFL, LASA laboratory
%   # Email: jrmout@gmail.com

d=size(data,1)/2;
options_solver=sdpsettings('solver',options.solver, ...
                           'verbose', options.verbose, ...
                           'warning', options.warning);

% Solver variables
A_inv = sdpvar(d,d,'full');
x_star = sdpvar(d,1);
error = sdpvar(d,size(data,2));
objective_function=sum((sum(error.^2)));
if isfield(options, 'weights')
    objective_function=sum(options.weights.*(sum(error.^2)));
end

% Define constraints
C=[error == -A_inv*data(d+1:2*d,:) + repmat(x_star,1,size(data,2))-data(1:d,:) ];
% Positive definite
C = C + [A_inv + A_inv' >= options.eps_constraints*eye(d,d)] ;

% Do not estimate the bias, set it to the one specified a priori
if options.bias == false
    if (size(options.attractor,1) ~= d)
        error(['The specified attractor should have size ' d 'x1']);
    else
        C = C + [x_star == options.attractor];
    end
end

% Solve the optimization
sol = optimize(C,objective_function,options_solver);
if sol.problem~=0
    warning(sol.info);
end
A_inv = value(A_inv);
x_attractor = value(x_star);
