% Setup paths
setup_stable_lds;

% Get trajectories from mouse
limits = [0 100 0 100];
data = generate_mouse_data(limits);
n_comp = 7;
em_iterations = 20;

% Optimization options
clear options;
options.n_iter = em_iterations;        % Max number of EM iterations
options.solver = 'fminsdp';              % Solver
options.criterion = 'mse';              % Solver
options.c_reg = -1e-4;                  % Pos def eps margin
options.verbose = 1;                    % Verbose (0-5)
options.warning = true;                % Display warning information
options.max_iter = 30;

lambda = em_mix_inv_lds(data, n_comp, options);

% Plot result
plot_streamlines_mix_lds_inv(lambda,limits);