function [cost, dcost_dp] = weighted_mse_mix(p, d, n_comp, data, weights)
% This function computes the weighted mean squared error of a linear system. 
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
    cost = cost + sum(weights(i,:).*(sum(error(:,:,i).^2)));
    if nargout > 1
        dcost_dA(:,:,i) = 2*(repmat(weights(i,:), [d 1]).*error(:,:,i))*data(1:d,:)';
        dcost_db = dcost_db + sum(2*(repmat(weights(i,:), [d 1]).*error(:,:,i)),2);
    end
end
if nargout > 1
    dcost_dp = fold_mix_lds(dcost_dA, dcost_db);
end