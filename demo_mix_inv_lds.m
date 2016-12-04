% Setup paths
setup_stable_lds;

% Get trajectories from mouse
limits = [0 100 0 100];
data = generate_mouse_data(limits);
n_comp = 5;
em_iterations = 5;

% Optimization options
clear options;
options.n_iter = 5;                     % Max number of EM iterations
options.solver = 'sedumi';              % Solver
options.eps_pos_def = 0.001;            % Pos def eps margin
options.verbose = 0;                    % Verbose (0-5)
options.warning = false;                % Display warning information

lambda = em_mix_inv_lds(data, n_comp, options);

% Plot result
plot_streamlines_mix_lds_inv(lambda,limits);