function [A_out, b_out]=estimate_stable_lds(data, varargin)
% ESTIMATE_STABLE_LDS fits a stable linear dynamical system x_dot = A*x + b
%   to data
%
%   ESTIMATE_STABLE_LDS solves the optimization problem
%    min (error^2)  subject to:  error == A*x + b - x_dot
%     A                          A'+A <= -options.eps_constraints*I
%   
%   [A_out, b_out] = ESTIMATE_STABLE_LDS(data, options) returns the system
%   matrix and the bias of a linear dynamical system
%
%   USAGE:
%   [A_out, b_out] = ESTIMATE_STABLE_LDS(data) fits a linear
%   dynamical system to the data and returns the system matrix
%   and the estimated bias.
%
%   [A_out, b_out] = ESTIMATE_STABLE_LDS(data, options) fits a 
%   linear dynamical system to the data with the specified options 
%   
%   INPUT PARAMETERS:
%   -data    data = [x; x_dot] and size(data) = [d*2,n_data_points], 
%            where d is the dimenstion of the input/output.
%   -options options.solver -- specifies the YALMIP solver or fmincon. 
%            options.eps_pos_def     -- specifies the eps for the
%                                       constraints
%            options.attractor       -- specifies a priori the atractor
%            options.weights         -- weighting factor for each sample
%            options.verbose         -- verbose YALMIP option [0-5]
%            options.warning         -- warning YALMIP option (true/false)
%
%   OUTPUT PARAMETERS:
%   - A_out  estimated system matrix
%   - b_out  estimated bias
%
%
%   # Authors: Jose Medina
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
if ~isfield(options, 'weights')
    options.weights =  ones(1,size(data,2));
end

d=size(data,1)/2;

if ~strcmp(options.solver, 'fmincon')
    %% YALMIP solvers
    options_solver=sdpsettings('solver',options.solver, ...
                               'verbose', options.verbose, ...
                               'warning', options.warning);
                           
    % Solver variables
    A = sdpvar(d,d,'full');
    b = sdpvar(d,1);
    error = sdpvar(d,size(data,2));

    % Objective
    objective_function=sum(options.weights.*(sum(error.^2)));

    % Constraints
    C=[error == (A*data(1:d,:) + repmat(b,1,size(data,2)))-data(d+1:2*d,:)];
    % Lyapunov LMI setting P=I -> Pos def (nonsymmetric) matrix
    C = C + [A'+A <= -options.eps_pos_def*eye(d,d)] ;

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
else
    %% MATLAB's fmincon
    opt_options = optimset( 'Algorithm', 'interior-point', 'LargeScale', 'off',...
        'GradObj', 'on', 'GradConstr', 'on', 'DerivativeCheck', 'off', ...
        'Display', 'iter', 'TolX', 1e-12, 'TolFun', 1e-12, 'TolCon', 1e-12, ...
        'MaxFunEval', 200000, 'MaxIter', 1000, 'DiffMinChange', ...
         1e-10, 'Hessian','off');
     
    A0 = eye(d);
    b0 = zeros(d,1);
    p0 = fold_lds(A0,b0);
    
    constraints_handle = @(p)neg_def_lmi(p, d, options);
    objective_handle = @(p)weighted_mse_linear(p, d, data, options.weights);

    % Solve
    p_opt = fmincon(objective_handle, p0, [], [], [], [], [], [], ...
                                              constraints_handle, opt_options);
    
	% Reshape A and b
    [A_out,b_out] = unfold_lds(p_opt,d);
end
