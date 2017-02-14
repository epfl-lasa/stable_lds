function [lambda]=em_mix_inv_lds(data, n_comp, varargin)
% EM_MIX_INV_LDS fits a Gaussian mixture of stable inverse linear 
% dynamical systems to data with the EM algorithm.
% The model assumed is given by
%  p([x_loc ; x_reg]) =  [0 ; x_attractor] + sum_i^n_comp pi_i * 
%                 N([mu_xloc ; -A^{-1}*x_dot], [Sigma_xloc 0 ; 0 Sigma_reg])
%   USAGE:
%   [lambda] = EM_MIX_INV_LDS(data, n_comp) fits a linear
%   dynamical system to the data and returns the inverse of the system matrix
%   and the estimated attractor.
%
%   [lambda] = EM_MIX_INV_LDS(data, n_comp) fits a linear
%   dynamical system to the data and returns the inverse of the system matrix
%   and the estimated attractor.
%
%   INPUT PARAMETERS:
%   -data    data = [x; x_dot] and size(data) = [d*2,n_data_points], 
%            where d is the dimenstion of the input/output.
%   -n_comp  number of components of the mixture.
%   -options options.n_iter -- number of iterations of the EM algorithm
%            options.min_eig_loc -- minimum eigenvalues for the covariance 
%                                   matrix of the local gaussians Sigma_xloc
%            options.min_eig_reg -- minimum eigenvalues for the covariance 
%                                   matrix of the regression noise Sigma_reg
%
%            Note that the same options variable is passed to the 
%            estimate_stable_mix_inv_lds function to be able to configure
%            any desired optimization parameters.
%
%   OUTPUT PARAMETERS:
%   - lambda    is a structure with the maximum likelihood parameters with
%               fields:
%                   - lambda.mu_xloc{n_comp} and lambda.cov_xloc{n_comp}
%                   that define the n_comp Gaussians of the x_loc variable
%                   - lambda.A_inv{n_comp}, lambda.cov_reg{n_comp} that
%                   define the n_comp inverse linear dynamical systems
%                   sharing the same attractor lambda.x_attractor
%                   - lambda.pi(n_comp) which represents the mixture weights 
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
if ~isfield(options, 'n_iter')
    options.n_iter = 5;
end
if ~isfield(options, 'min_eig_loc')
    options.min_eig_loc = 1e-100;
end
if ~isfield(options, 'min_eig_reg')
    options.min_eig_reg = 1e-100;
end

d = size(data,1)/2;
x_obs = data(1:d,:)';
x_dot_obs = data(d+1:end,:)';

% Init model parameters with kmeans
lambda= init_kmeans_mix_lds(data, n_comp, options);
weights = zeros(n_comp, size(data,2));
loglik = zeros(options.n_iter,1);

for it = 1:options.n_iter
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
    loglik(it) = sum(log(sum(weights,1)));
    weights = weights ./ repmat(sum(weights,1), n_comp, 1);
    
    if it > 1 && (abs(loglik(it) - loglik(it-1)) < 1e-3) 
        break;
    end

    %% M step
    % Initial values for the optimization
    for i = 1:n_comp
        A_inv_old(:,:,i) = -lambda.A_inv{i};
    end
    % Max A_invs and x_attractor    
    [lambda.x_attractor, lambda.A_inv] = estimate_stable_mix_inv_lds( ...
                                [x_obs x_dot_obs]', weights, options, lambda.x_attractor, A_inv_old);

    for c=1:n_comp
        % Max regression error covariance cov_reg
        model_error = (-lambda.A_inv{c}*x_dot_obs' ...
                 + repmat(lambda.x_attractor, [1 size(x_obs,1)]) - x_obs');
        cov_reg = 1/sum(weights(c,:))* ...
                        (repmat(weights(c,:), size(model_error,1), 1) ... 
                            .* model_error*model_error');
        lambda.cov_reg{c} = crop_min_eig(cov_reg, options.min_eig_reg);

        % Max local gaussian x_loc
        w_factor = weights(c,:);
        mu_c_loc = 1/(sum(w_factor)) * ...
                              sum((repmat(w_factor,size(x_obs,2),1)'.*x_obs));

        % Covariance
        dev = (repmat(mu_c_loc, size(x_obs,1), 1) - x_obs);
        cov_c_loc =  1/(sum(w_factor)) * ...
                                 dev'*(repmat(w_factor,size(dev,2),1)'.*dev);
                             
        lambda.cov_xloc{c} = crop_min_eig(cov_c_loc, options.min_eig_loc);
        lambda.mu_xloc{c} = mu_c_loc';
    end
    lambda.pi = sum(weights,2)/sum(sum(weights));
end
loglik
end