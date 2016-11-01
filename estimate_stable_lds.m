function [A_out, b_out]=estimate_stable_lds(data, options, varargin)
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
%   -varargin{1} specifies a priori the attractor of the DS
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
options_solver=sdpsettings('solver',options.solver, ...
                           'verbose', options.verbose, ...
                           'warning', options.warning);

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
% Set the attractor at the specified input if set a priori
if nargin > 2
    % Do not estimate the bias, set it to the one specified a priori
    if options.bias == false
        if (size(varargin{1},1) ~= d)
            error(['The specified attractor should have size ' d 'x1']);
        else
            C = C + [A*varargin{1}-b == zeros(d,1)];
        end
    end
end

% Solve the optimization
sol = optimize(C,objective_function,options_solver);
if sol.problem~=0
    warning(sol.info);
end

A_out = value(A);
b_out = value(b);
