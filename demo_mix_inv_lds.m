% Setup paths
setup_stable_lds;

% Get trajectories from mouse
limits = [0 100 0 100];
data = generate_mouse_data(limits);
n_comp = 3;
em_iterations = 5;

% Optimization options
options.solver = 'mosek';
options.eps_pos_def = 0.1;
options.verbose = 1; 
options.warning = 0;

lambda = em_mix_lds_inverse(data, n_comp, em_iterations, options);

% Plot result
plot_streamlines_mix_lds_inv(lambda,limits);