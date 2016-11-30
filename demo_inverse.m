% Setup paths
setup_stable_lds;

% Get trajectories from mouse
limits = [0 100 0 100];
data = generate_mouse_data(limits);

options.solver = 'mosek';
options.eps_constraints = 0.01;
options.bias = true;
options.verbose = 1;
options.warning = false;
options.attractor = [0 0]';
options.weights = ones(1,size(data,2)); % normalized weights for each sample

[A_inv, x_attractor] = estimate_stable_lds_inverse(data,options);

% Plot result
plot_streamlines_inverse(A_inv, x_attractor, limits);