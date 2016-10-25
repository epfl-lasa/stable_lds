% Setup paths
setup_stable_lds;

% Get trajectories from mouse
data = generate_mouse_data();

options.solver = 'sedumi';
options.lyapunov_lmi = 'simplified';
options.eps_constraints = 1e-4;
options.bias = true;

[A, b, P] = estimate_stable_lds(data,options);

% Plot result
