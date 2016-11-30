function [ lambda ] = init_kmeans_mix_lds( data, n_comp, options)
%INIT_KMEANS_MIX_LDS Initializes the model with kmeans

d=size(data,1)/2;
x_obs = data(1:d,:)';
x_dot_obs = data(d+1:end,:)';

% The structure containing all model parameters
lambda.mu_xloc = cell(n_comp,1);
lambda.cov_xloc = cell(n_comp,1);
lambda.cov_reg = cell(n_comp,1);

% Compute kmeans
[idx,~] = kmeans([x_obs x_dot_obs], n_comp, 'Replicates',10);

% Fit a stable LDS to each of the kmeans clusters
options.solver = 'sedumi';
options.eps_pos_def = 0.01;
options.verbose = 1; 
options.warning = 0;

weights = zeros(n_comp, size(data,2));
for c=1:n_comp
    weights(c,:) = (idx == c)';
end

[lambda.x_attractor, lambda.A_inv] = estimate_stable_mix_inv_lds( ...
                                 [x_obs x_dot_obs]', weights, n_comp, options);

for c=1:n_comp
    x_obs_c = x_obs(idx == c,:);
    x_dot_obs_c = x_dot_obs(idx == c,:);

    % Estimate noise from prediction error covariance
    model_error = (-lambda.A_inv{c}*x_dot_obs_c' ...
      + repmat(lambda.x_attractor, [1 size(x_obs_c,1)]) - x_obs_c');
    lambda.cov_reg{c} = 1/size(x_obs_c,2)*(model_error*model_error');
  
    % Compute the variance of the inputs
    % lambda.cov_reg{c} = diag(diag(cov(x_obs_c)));
    
    lambda.mu_xloc{c} = mean(x_obs_c)';
    lambda.cov_xloc{c} = cov(x_obs_c);
end
lambda.pi = (1/n_comp) * ones(n_comp,1);

end

