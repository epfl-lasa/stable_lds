function x_dot = get_dyn_mix_inv_lds(lambda, x)
n_comp = length(lambda.pi);
x_dot = zeros(size(x));

weights = zeros(n_comp, size(x,2));

% Compute dynamics
for c=1:n_comp
    weights(c,:) = ( mvnpdf(x', lambda.mu_xloc{c}', ...
                                         lambda.cov_xloc{c}) ...
                   .* lambda.pi(c) )' + realmin; % TODO: In case of numerical 
                                                 % problems, instead of realmin
                                                 % choose closest component 
                                                 % based on Mahalanobis dist
end
weights = weights ./ (repmat(sum(weights,1), n_comp, 1) + n_comp*realmin); 

sum_A_inv = zeros(2,2,length(x));
for c=1:n_comp
    if (sum( eig(lambda.A_inv{c} + lambda.A_inv{c}') <= 0 ) > 0)
        eig(lambda.A_inv{c} + lambda.A_inv{c}')
        disp('Not stable!')
    end
    sum_A_inv = sum_A_inv + repmat(reshape(weights(c,:), ...
                                [1 1 length(weights(c,:))]), 2,2,1)...
                                .*repmat(lambda.A_inv{c},1,1,length(x));
end
for i=1:length(x)
    x_dot(:,i) = -sum_A_inv(:,:,i)\(x(:,i) - lambda.x_attractor);
end