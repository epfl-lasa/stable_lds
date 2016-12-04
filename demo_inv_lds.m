% Setup paths
setup_stable_lds;

% Get trajectories from mouse
limits = [0 100 0 100];
disp('Draw a few trajectories in the figure with your mouse and press stop recording when done.')
data = generate_mouse_data(limits);

% Optimization options
clear options;
options.solver = 'sedumi';
options.eps_constraints = 0.001;        % Pos def eps margin
options.verbose = 0;                    % Verbose (0-5)
options.warning = false;                % Display warning information
%options.attractor = [0 0]';             % Ignored if options.bias=true
options.weights = ones(1,size(data,2)); % Weights for each sample

[A_inv, x_attractor] = estimate_stable_inv_lds(data, options);

% Plot result
plot_streamlines_inv_lds(A_inv, x_attractor, limits);