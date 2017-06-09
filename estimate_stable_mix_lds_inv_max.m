function [x_attractor, A_out]=estimate_stable_mix_lds_inv_max(data, ...
                                                            weights, varargin)
% ESTIMATE_STABLE_MIX_LDS_INV_MAX fits weighted sum of n_comp stable linear
% dynamical systems to a weighted sum of datapoints. The optimization
% problem is for mse criterion 
%
% min sum_{c=1}^{n_comp} (weights(c,:).*||(x - (x_star - A_inv_c*x_dot))||
%
% or for the maximum likelihood criterion (logdet)
%
% min sum_{c=1}^{n_comp} log( det ( weights(c,:).*
%         (x - (x_star - A_inv_c*x_dot)) * (x - (x_star -A_inv_c*x_dot))' ) )
%
% The logdet criterion is necessary for maximum likelihood estimation
% methods.
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
%
%   Optional input parameters
%   -options options.solver -- specifies the YALMIP solver or fmincon. 
%            options.criterion       -- mse|logdet      
%            options.weights         -- weighting factor for each sample
%            options.verbose         -- verbose YALMIP option [0-5]
%
%            %% NLP solvers only
%            options.max_iter        -- Max number of iterations
%            options.c_reg           -- specifies the eps constant for the
%                                       constraint
%
%            %% YALMIP solvers only
%            options.eps_pos_def     -- specifies the eps constant for the
%                                       constraints
%            options.warning         -- warning YALMIP option (true/false)
%            options.attractor       -- specifies a priori the atractor
%
%   This code provides 3 different solvers for this problem
%   - YALMIP: solvers for convex problems
%   - fmincon: NLP solver with a nonconvex constraint for the LMI
%   - fminsdp: NLP solver with a convex constraint for the LMI solved with
%              the ldl method
%
%   In general any YALMIP solver will provide the best solution and should 
%   be the preferred option. 
%   The NLP solver works fine most of the time but may provide 
%   sometimes nonstable solutions (can't find feasible ones).
%   The NLP solutions are there only as a first step to combine this
%   convex constrained problem with nonlinear nonconvex settings without
%   adding any nonconvexity to the problem.
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
    if nargin > 3
        A_0 = varargin{3};
        x_star_0 = varargin{2};
    end
else
    options = [];
end

% Default values
if ~isfield(options, 'solver')
    options.solver = 'fminsdp';
end 
if ~isfield(options, 'verbose')
    options.verbose = 0;
end 
if ~isfield(options, 'criterion')
    options.criterion = 'logdet';
end 

n_comp = size(weights,1);
A_out = cell(n_comp,1);
d=size(data,1)/2;

%% YALMIP SQP solvers
if ~isfield(options, 'c_reg')
    options.c_reg = 1e-3;
end
if ~isfield(options, 'warning')
    options.warning = 0;
end
if ~isfield(options, 'debug')
    options.debug = 0;
end

options_solver=sdpsettings('solver',options.solver, ...
                           'verbose', options.verbose, ...
                           'warning', options.warning, ...
                           'debug', options.debug);

% Do not estimate the attractor, set it to the one specified a priori
if ~isfield(options, 'attractor')
    yalmip('clear');
    % Solver variables
    A_inv = sdpvar(d,d,n_comp,'full');
    x_star = sdpvar(d,1);
    error = sdpvar(d,size(data,2), n_comp);
    objective_function=0;
    C = [];
    for i = 1:n_comp
        objective_function = objective_function ...
         + (1/size(data,2))*sum(weights(i,:))*((1/size(data,2)) * ...
                                    sum(weights(i,:).*(sum(error(:,:,i).^2))));
        C = C + [error(:,:,i) == -A_inv(:,:,i)*data(d+1:2*d,:) ...
                                 + repmat(x_star,1,size(data,2))-data(1:d,:) ];
        C = C + [A_inv(:,:,i) + A_inv(:,:,i)' >= options.c_reg_inv*eye(d,d)];
    end
    if isfield(options, 'prior')
        objective_function = objective_function + ...
            0.5*(x_star - options.prior.mu)' * options.prior.sigma_inv * ...
                                        (x_star - options.prior.mu);
    end

    % Solve the optimization
    sol = optimize(C,objective_function,options_solver);
    if sol.problem~=0
        warning(sol.info);
    end
    x_attractor = value(x_star);

else
    x_attractor = options.attractor;
end
yalmip('clear');
% Solver variables
A = sdpvar(d,d,n_comp,'full');
error = sdpvar(d,size(data,2), n_comp);
objective_function=0;
C = [];

for i = 1:n_comp
    objective_function = objective_function ...
     + (1/size(data,2))*sum(weights(i,:))*((1/size(data,2)) * ...
                                sum(weights(i,:).*(sum(error(:,:,i).^2))));
    C = C + [error(:,:,i) == A(:,:,i)*(data(1:d,:) ...
                   - repmat(x_attractor,1,size(data,2)))-data(d+1:2*d,:) ];
    C = C + [A(:,:,i) + A(:,:,i)' <= -options.c_reg*eye(d,d)];
end
% Solve the optimization
sol = optimize(C,objective_function,options_solver);
if sol.problem~=0
    warning(sol.info);
end

for i = 1:n_comp
    A_out{i} = value(A(:,:,i));
end


end
