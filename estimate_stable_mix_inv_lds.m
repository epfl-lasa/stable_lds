function [x_attractor, A_inv_out]=estimate_stable_mix_inv_lds(data, ...
                                                            weights, varargin)
% ESTIMATE_STABLE_MIX_INV_LDS fits weighted sum of n_comp stable inverse linear
% dynamical systems to a weighted sum of datapoints. The optimization
% problem is 
%
% min sum_{c=1}^{n_comp} ||(weights(c,:).*(x - (x_star - A_inv_c*x_dot)).^2)||
%
%   USAGE:
%   [A_inv, x_attractor] = ESTIMATE_STABLE_MIX_INV_LDS(data, weights) fits 
%    a mixture of inverse linear dynamical system to the data with he
%    corresponding weights and returns the inverse system matrices and the
%    attractor
%
%   [A_inv, x_attractor] = ESTIMATE_STABLE_MIX_INV_LDS(data, weights, options) 
%   considers also the optional parameters from options
% 
%   INPUT PARAMETERS:
%   -data    data = [x; x_dot] and size(data) = [d*2,n_data_points], 
%            where d is the dimenstion of the input/output.
%   -weights weights for each component for each datapoint
%            size(weights) = [n_comp, n_data_points]
%   -options options.solver          -- specifies the YALMIP solver.
%            options.eps_pos_def     -- specifies the eps for the
%                                       positive definiteness constraint
%            options.attractor       -- specifies the attractor a priori
%            options.verbose         -- verbose YALMIP option [0-5]
%            options.warning         -- warning YALMIP option (true/false)
%
%   OUTPUT PARAMETERS:
%   - x_attractor  estimated attractor of the system
%   - A_inv_out    cell(n_comp x 1) with the inverse system matrices of the
%                  model
%
%   # Author: Jose Medina
%   # EPFL, LASA laboratory
%   # Email: jrmout@gmail.com

% Check for options
if nargin > 2
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
if ~isfield(options, 'eps_pos_def')
    options.eps_pos_def = 1e-3;
end

n_comp = size(weights,1);
A_inv_out = cell(n_comp,1);
d=size(data,1)/2;

if ~strcmp(options.solver, 'fmincon')
    %% YALMIP solvers
    options_solver=sdpsettings('solver',options.solver, ...
                               'verbose', options.verbose, ...
                               'warning', options.warning);

    % Solver variables
    A_inv = sdpvar(d,d,n_comp,'full');
    x_star = sdpvar(d,1);
    error = sdpvar(d,size(data,2), n_comp);
    objective_function=0;
    C = [];
    for i = 1:n_comp
        objective_function = objective_function ...
                                    + sum(weights(i,:).*(sum(error(:,:,i).^2)));
        C = C + [error(:,:,i) == -A_inv(:,:,i)*data(d+1:2*d,:) ...
                                    + repmat(x_star,1,size(data,2))-data(1:d,:) ];
        C = C + [A_inv(:,:,i) + A_inv(:,:,i)' >= options.eps_pos_def*eye(d,d)];
    end

    % Do not estimate the attractor, set it to the one specified a priori
    if isfield(options, 'attractor')
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

    for i = 1:n_comp
        A_inv_out{i} = value(A_inv(:,:,i));
    end
    x_attractor = value(x_star);
else
     %% MATLAB's fmincon
    opt_options = optimset( 'Algorithm', 'interior-point', 'LargeScale', 'off',...
        'GradObj', 'on', 'GradConstr', 'on', 'DerivativeCheck', 'off', ...
        'Display', 'iter', 'TolX', 1e-20, 'TolFun', 1e-20, 'TolCon', 1e-10, ...
        'MaxFunEval', 200000, 'MaxIter', 1000, 'DiffMinChange', ...
         1e-10, 'Hessian','off');
    
    x_star_0 = mean(data(1:d,:)')';
    for i = 1:n_comp
        A_inv_0(:,:,i)= eye(d);
    end
    p0 = fold_mix_lds(A_inv_0,x_star_0);
    
    data_inv = zeros(size(data));
    data_inv(1:d,:) = data(d+1:end,:);
    data_inv(d+1:end,:) = data(1:d,:);
    
    constraints_handle = @(p)neg_def_lmi_mix(p, d, n_comp, options);
    objective_handle = @(p)weighted_mse_mix(p, d, n_comp, data_inv, weights);

    % Solve
    p_opt = fmincon(objective_handle, p0, [], [], [], [], [], [], ...
                                              constraints_handle, opt_options);
	% Reshape A and b
    [A_inv_neg,x_attractor] = unfold_mix_lds(p_opt,d,n_comp);
    for i = 1:n_comp
        A_inv_out{i} = -A_inv_neg(:,:,i);
    end
end