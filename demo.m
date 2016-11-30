% Setup paths
close all;
setup_stable_lds;

% Get trajectories from mouse
limits = [0 100 0 100];
data = generate_mouse_data(limits);

options.solver = 'mosek';
options.eps_constraints = 0.1;
options.bias = true;
options.verbose = 5;
options.warning = true;
options.criterion = 'mse';
options.attractor = [0 0]';
options.weights = ones(1,size(data,2)); % normalized weights for each sample


figure(2);
options.eps_constraints = 0;
[A_direct, b_direct] = estimate_stable_lds(data,options);
plot_streamlines(A_direct, b_direct, limits);