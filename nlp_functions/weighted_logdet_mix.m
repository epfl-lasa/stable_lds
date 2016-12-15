function [cost] = weighted_logdet_mix(p, d, n_comp, data, weights)
% This function computes the weighted logdet covariance of a mixture of linear 
% systems. 
% It returns the value and the derivatives w.r.t. the bias and the
% linear term
[A,b] = unfold_mix_lds(p,d,n_comp);
cost=0;
error = zeros(d,size(data,2),n_comp);
dcost_db = zeros(d,1);
dcost_dA = zeros(d,d,n_comp);

for i = 1:n_comp
    error(:,:,i) = A(:,:,i)*data(1:d,:) ...
                                + repmat(b,1,size(data,2))-data(d+1:2*d,:); 
    covariance = 1/sum(weights(i,:))* ...
        (repmat(weights(i,:), d, 1).*error(:,:,i)*error(:,:,i)');
    cost = cost + log(det(covariance));
end
%% TODO: Add derivatives. Look for an efficient way of computing it.