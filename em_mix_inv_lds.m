function [lambda]=em_mix_inv_lds(data, n_comp, n_iter, options)
% EM_MIX_LDS_INVERSE fits a Gaussian mixture of stable inverse linear 
% dynamical systems to data with the EM algorithm.
% The model assumed is given by
%  p([x_loc ; x_reg]) =  [0 ; x_attractor] + sum_i^n_comp pi_i * 
%                 N([mu_xloc ; -A^{-1}*x_dot], [Sigma_xloc 0 ; 0 Sigma_reg])
%
%   INPUT PARAMETERS:
%   -data    data = [x; x_dot] and size(data) = [d*2,n_data_points], 
%            where d is the dimenstion of the input/output.
%   -n_comp  number of components of the mixture.
%   -options options.solver -- specifies the YALMIP solver
%            options.eps_constraints -- specifies the eps for the
%                                       pos def constraint satisfaction
%
%   OUTPUT PARAMETERS:
%   - model_params    is a structure with the maximum likelihood parameters
%
%   # Author: Jose Medina
%   # EPFL, LASA laboratory
%   # Email: jrmout@gmail.com

d = size(data,1)/2;
x_obs = data(1:d,:)';
x_dot_obs = data(d+1:end,:)';

min_eig_loc = 1;
min_eig_reg = 1e-1;

% Init model parameters with kmeans
lambda= init_kmeans_mix_lds(data, n_comp, min_eig_reg, min_eig_loc, options);

weights = zeros(n_comp, size(data,2));

for it = 1:n_iter
    %% E step    
    for c=1:n_comp
        weights(c,:) = ( mvnpdf(x_obs,...
                        (repmat(lambda.x_attractor, 1, size(data,2)) ...
                                        - lambda.A_inv{c}*x_dot_obs')',...
                         lambda.cov_reg{c}) ...
                       .* mvnpdf(x_obs, lambda.mu_xloc{c}', ...
                                             lambda.cov_xloc{c}) ...
                       .* lambda.pi(c) )';
    end
    weights = weights ./ repmat(sum(weights,1) + 1e-10, n_comp, 1);

    %% M step
    % Max A_invs and x_attractor
    [lambda.x_attractor, lambda.A_inv] = estimate_stable_mix_inv_lds( ...
                                [x_obs x_dot_obs]', weights, n_comp, options);

    for c=1:n_comp
        % Max regression error covariance cov_reg
        model_error = (-lambda.A_inv{c}*x_dot_obs' ...
                 + repmat(lambda.x_attractor, [1 size(x_obs,1)]) - x_obs');
        cov_reg = 1/sum(weights(c,:))* ...
                        (repmat(weights(c,:), size(model_error,1), 1) ... 
                            .* model_error*model_error');
        lambda.cov_reg{c} = diag(diag(crop_min_eig(cov_reg, min_eig_reg)));

        % Max local gaussian x_loc
        w_factor = weights(c,:);
        mu_c_loc = 1/(sum(w_factor)) * ...
                              sum((repmat(w_factor,size(x_obs,2),1)'.*x_obs));

        % Covariance
        dev = (repmat(mu_c_loc, size(x_obs,1), 1) - x_obs);
        cov_c_loc =  1/(sum(w_factor)) * ...
                                 dev'*(repmat(w_factor,size(dev,2),1)'.*dev);
                             
        lambda.cov_xloc{c} = crop_min_eig(cov_c_loc, min_eig_loc);
        lambda.mu_xloc{c} = mu_c_loc';
    end
    lambda.pi = (1/n_comp) * ones(n_comp,1);
end
end