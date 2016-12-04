function [A_inv, x_attractor]=estimate_stable_inv_lds(data, varargin)
%   ESTIMATE_STABLE_INV_LDS fits a stable linear dynamical system inverting
%   the estimation problem. Instead of considering the standard model
%       x_dot = A*x + b = -A (x - x_attractor)
%   it assumes that A is positive definite (and possibly nonsymmetric) and 
%   performs the regression task in the inverse model, i.e.
%       x = x_attractor - A_inv x_dot
%   this model has the benefit of having the attractor as bias. More
%   precisely it solves the problem
%   min (error^2)  subject to:  error == x_attractor - A_inv x_dot - x
%       A                       A_inv + A_inv' >= options.eps_constraints*I
%
%   USAGE:
%   [A_inv, x_attractor] = ESTIMATE_STABLE_INV_LDS(data) fits a linear
%   dynamical system to the data and returns the inverse of the system matrix
%   and the estimated attractor.
%
%   [A_inv, x_attractor] = ESTIMATE_STABLE_INV_LDS(data, options) fits a 
%   linear dynamical system to the data with the specified options 
%   
%   INPUT PARAMETERS:
%   -data    data = [x; x_dot] and size(data) = [d*2,n_data_points], 
%            where d is the dimenstion of the input/output.
%
%   Optional input parameters
%   -options options.solver -- specifies the YALMIP solver.
%            options.eps_constraints -- specifies the eps for the
%                                       constraints
%            options.attractor       -- specifies a priori the atractor
%            options.weights         -- weighting factor for each sample
%
%   OUTPUT PARAMETERS:
%   - A_inv  estimated inverse system matrix
%   - x_attractor  estimated attractor
%
%   # Author: Jose Medina
%   # EPFL, LASA laboratory
%   # Email: jrmout@gmail.com

% Check for options
if nargin > 1
    options = varargin{1};
else
    options = [];
end

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
A_inv = sdpvar(d,d,'full');
x_star = sdpvar(d,1);
error = sdpvar(d,size(data,2));
objective_function=sum((sum(error.^2)));
if isfield(options, 'weights')
    objective_function=sum(options.weights.*(sum(error.^2)));
end

% Constraints
C=[error == -A_inv*data(d+1:2*d,:) + repmat(x_star,1,size(data,2))-data(1:d,:) ];
% Lyapunov LMI setting P=I -> Pos def (nonsymmetric) matrix
C = C + [A_inv + A_inv' >= options.eps_constraints*eye(d,d)] ;

% Do not estimate the bias, set it to the one specified a priori
if isfield(options, 'attractor')
    if (size(options.attractor,1) ~= d)
        error(['The specified attractor should have size ' d 'x1']);
    else
        C = C + [x_star == options.attractor];
    end
end

% Solve
sol = optimize(C,objective_function,options_solver);
if sol.problem~=0
    warning(sol.info);
end
A_inv = value(A_inv);
x_attractor = value(x_star);
