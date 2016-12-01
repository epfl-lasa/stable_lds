function [x_attractor, A_inv_out]=estimate_stable_mix_inv_lds(data, ...
                                                    weights, n_comp, options)
% ESTIMATE_STABLE_MIX_INV_LDS fits weighted sum of stable linear
% dynamical systems to a weighted sum of datapoints. The optimization
% problem is 
%
% min sum_{c=1}^{n_comp} ||(weights(:,c).*(x - (x_star - A_inv_c*x_dot)))||
%
%   
%   INPUT PARAMETERS:
%   -data    data = [x; x_dot] and size(data) = [d*2,n_data_points], 
%            where d is the dimenstion of the input/output.
%   -weights weights for each component for each datapoint
%            size(weights) = [n_comp, n_data_points]
%   -options options.solver -- specifies the YALMIP solver.
%            options.eps_pos_def -- specifies the eps for the
%                                   positive definiteness constraint
%            options.attractor -- (optional) specifies the attractor a priori
%
%   OUTPUT PARAMETERS:
%   - lambda  structure with the estimated parameters
%
%   # Author: Jose Medina
%   # EPFL, LASA laboratory
%   # Email: jrmout@gmail.com

d=size(data,1)/2;
options_solver=sdpsettings('solver',options.solver, ...
                           'verbose', options.verbose, ...
                           'warning', options.warning);
A_inv_out = cell(n_comp,1);

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
    eig(A_inv_out{i})
end
x_attractor = value(x_star);