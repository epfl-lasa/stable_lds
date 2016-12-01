function cov_cropped = crop_min_eig(cov, min_eig)
% CROP_MIN_EIG(COV,MIN_EIG) crops covariance eigenvalues to a minimum value
    [V,D] = eig(cov);
    d_eig = diag(D);
    d_eig(d_eig<min_eig) = min_eig;
    cov_cropped = V*diag(d_eig)*V';
end

