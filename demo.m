% Setup paths
setup_stable_lds;

% Get trajectories from mouse
limits = [0 100 0 100];
data = generate_mouse_data(limits);

options.solver = 'sedumi';
options.eps_constraints = 1e-4;
options.bias = true;

[A, b] = estimate_stable_lds(data,options);

% Plot result
plot_streamlines(A, b, limits)