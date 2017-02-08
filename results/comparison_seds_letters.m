% NOTE: To run this comparison, clone the seds repository first in the
% source directory of this repo. https://bitbucket.org/khansari/seds

run ../setup_stable_lds;

% For Mosek solver (use sedumi if you don't have a license)
% addpath('/home/medina/Dropbox/work/3rdParty/mosek/8/toolbox/r2014a')

condition = '3';
do_plot = false;
max_c = 20; % Maximum number of components to evaluate 1:max_c


switch condition
    case '1'
        seds_objective = 'mse';
        filename = 'results/results_fixed_attractor_mse';
        attractor_fixed = true;
    case '2'
        seds_objective = 'mse';
        filename = 'results/results_free_attractor_mse';
        attractor_fixed = false;
    case '3'
        seds_objective = 'likelihood';
        filename = 'results/results_fixed_attractor_likelihood';
        attractor_fixed = true;
end

files = {'Angle', 'Bump', 'CShape', 'GShape', 'JShape_2', 'JShape', ...
    'Khamesh', 'Line', 'Multi_Models_1', 'Multi_Models_2', 'Multi_Models_3', ...
    'Multi_Models_4', 'NShape', 'PShape', 'Rshape', 'Saeghe', 'Sharpc', ...
    'Sine', 'Soft_Sine', 'Spoon', 'Sshape', 'Trapezoid', 'WShape', 'Zshape'};

%% Optimization options SEDS
clear options;
dt = 0.1;
tol_cutting = 1;
% A set of options that will be passed to the solver. Please type 
% 'doc preprocess_demos' in the MATLAB command window to get detailed
% information about other possible options.
options_seds.tol_mat_bias = 10^-6; % A very small positive scalar to avoid
                              % instabilities in Gaussian kernel [default: 10^-15]
                              
options_seds.display = 1;          % An option to control whether the algorithm
                              % displays the output of each iterations [default: true]
                              
options_seds.tol_stopping=10^-10;  % A small positive scalar defining the stoppping
                              % tolerance for the optimization solver [default: 10^-10]

options_seds.max_iter = 1000;       % Maximum number of iteration for the solver [default: i_max=1000]

options_seds.objective = seds_objective;    % 'likelihood': use likelihood as criterion to
                              % optimize parameters of GMM
                              % 'mse': use mean square error as criterion to
                              % optimize parameters of GMM
                              % 'direction': minimize the angle between the
                              % estimations and demonstrations (the velocity part)
                              % to optimize parameters of GMM                              
                              % [default: 'mse']     

%% Optimization options SIEDS
clear options_sieds;
options_sieds.n_iter = 50;        % Max number of EM iterations
options_sieds.solver = 'sedumi';              % Solver 
options_sieds.criterion = 'mse';              % Solver
options_sieds.c_reg = 1e0;                  % Pos def eps margin
options_sieds.verbose = 3;                    % Verbose (0-5)
options_sieds.warning = true;                % Display warning information

if attractor_fixed
    options_sieds.attractor = [0;0];
end

mse_seds = zeros(max_c,1);
mse_sieds = zeros(max_c,1);

for c = 1:max_c
    for i=1:length(files)
    load(['models/recorded_motions/' files{i}],'demos');
    % the variable 'demos' composed of 3 demosntrations. Each demonstrations is
    % recorded from Tablet-PC at 50Hz. Datas are in millimeters.

    %% SEDS learning
    [x0 , xT, Data, index] = preprocess_demos(demos,dt,tol_cutting); %preprocessing datas
    [Priors_0, Mu_0, Sigma_0] = initialize_SEDS(Data,c); %finding an initial guess for GMM's parameter
    [Priors Mu Sigma]=SEDS_Solver(Priors_0,Mu_0,Sigma_0,Data,options_seds); %running SEDS optimization solver

    % A set of options that will be passed to the Simulator. Please type 
    % 'doc preprocess_demos' in the MATLAB command window to get detailed
    % information about each option.
    d = size(Data,1)/2; %dimension of data

    get_expected_dynamics_seds = @(x) GMR(Priors,Mu,Sigma,x,1:d,d+1:2*d);

    if do_plot
        % Plot
        figure('name','Streamlines','position',[800   90   560   320])
        plot(Data(1,:),Data(2,:),'r.')
        ax = gca;
        limits = [ax.XLim ax.YLim];
        plotStreamLines(Priors,Mu,Sigma,limits)
        hold on
        plot(Data(1,:),Data(2,:),'r.')
        title('Streamlines SEDS')
    end

    x_dot_seds = get_expected_dynamics_seds(Data(1:d,:));
    mse_seds(c) = mse_seds(c) + (1/size(Data,2)) * ...
                                    sum(sqrt(sum((x_dot_seds - Data(d+1:end, :)).^2)));

    %% SIEDS learning
    lambda = em_mix_inv_lds(Data, c, options_sieds);

    if do_plot
        figure;
        % Plot result
        plot(Data(1,:),Data(2,:),'r.');
        hold on;
        ax = gca;
        limits = [ax.XLim ax.YLim];
        plot_streamlines_mix_lds_inv(lambda,limits);
        title('Streamlines SIEDS')
    end

    x_dot_sieds = get_expected_dynamics(lambda, Data(1:d,:));
    mse_sieds(c) = mse_sieds(c) + (1/size(Data,2)) * ...
                                    sum(sqrt(sum((x_dot_sieds - Data(d+1:end, :)).^2)));
    end
    mse_seds(c) = mse_seds(c) / length(files);
    mse_sieds(c) = mse_sieds(c) / length(files);
end
save(filename, 'mse_seds', 'mse_sieds', 'options_seds', 'options_sieds');




