function run_demo_gui()
% A very simple self-explanatory GUI to quickly test the performance of
% the SIEDS model

setup_stable_lds;
% For mosek solver
addpath('~/Dropbox/work/3rdParty/mosek/8/toolbox/r2014a')

%% Params
% Model options
n_comp = 7;
em_iterations = 5;

% Optimization options
clear options;
options.n_iter = em_iterations;        % Max number of EM iterations
options.solver = 'mosek';              % Solver
options.criterion = 'mse';              % Solver
options.c_reg = 3e-1;                  % Pos def eps margin
options.verbose = 1;                    % Verbose (0-5)
options.warning = true;                % Display warning information
options.max_iter = 30;
limits = [0 100 0 100];
lambda = [];
p_handle = [];

% Window size for the Savitzky-Golay filter
f_window = 11;

%% Figure setup
fig = figure();
axes('Parent',fig,...
    'Position',[0.13 0.163846153846154 0.775 0.815]);
axis(limits);
hold on;
disp('Draw some trajectories with the mouse on the figure.')

% to store the data
X = [];         % unfiltered data
data = [];      % filtered data
dem_index = 0; % counter for demonstrations

% disable any figure modes
zoom off
rotate3d off
pan off
brush off
datacursormode off

% Setup callbacks for data capture
set(fig,'WindowButtonDownFcn',@(h,e)button_clicked(h,e));
set(fig,'WindowButtonUpFcn',[]);
set(fig,'WindowButtonMotionFcn',[]);
set(fig,'Pointer','circle');
hp = gobjects(0);

% Train button
uicontrol('style','pushbutton','String', 'Train','Callback',@train_recorded_motions, ...
          'position',[50 15 110 25], ...
          'UserData', 1);

% Clear button
uicontrol('style','pushbutton','String', 'Clear','Callback',@clear_trajectories, ...
          'position',[400 15 110 25], ...
          'UserData', 1);

      

%% Train button function
function train_recorded_motions(ObjectS, ~)
    set(ObjectS, 'UserData', 0); % unclick button
    % Train model
    lambda = em_mix_inv_lds(data, n_comp, options);
    delete(p_handle);
    [p_handle, l_handle] = plot_streamlines_mix_lds_inv(lambda,limits);
    set(l_handle, 'Position', [0.365 0.0191 0.291 0.124]);
end

%% Clear trajectories button function
function clear_trajectories(ObjectS, ~)
    set(ObjectS, 'UserData', 0); % unclick button
    delete(hp);
    delete(p_handle);
    X = [];
    data = [];
    dem_index = 0;
end

%% Functions for data capture
function ret = button_clicked(~,~)
    if(strcmp(get(gcf,'SelectionType'),'normal'));
        start_demonstration();
    end
end

function ret = start_demonstration()
    disp('Started demonstration');
    set(gcf,'WindowButtonUpFcn',@stop_demonstration);
    set(gcf,'WindowButtonMotionFcn',@record_current_point);
    ret = 1;
    tic;
end

function ret = stop_demonstration(~,~)
    disp('Stopped demonstration. Press train to learn a model with this data.');
    set(gcf,'WindowButtonMotionFcn',[]);
    set(gcf,'WindowButtonUpFcn',[]);
    set(gcf,'WindowButtonDownFcn',[]);
    dem_index = dem_index + 1;
    
    x_obs{dem_index} = X;
    X = [];
    
    % Savitzky-Golay filter and derivatives
    x_obs_dem = x_obs{dem_index}(1:2,:)';
    dt = mean(diff(x_obs{dem_index}(3,:)')); % Average sample ...
                                       % time (Third dimension contains time)
	if (size(x_obs_dem,1) > f_window)
        dx_nth = sgolay_time_derivatives(x_obs_dem, dt, 2, 3, f_window);
        data = [data [dx_nth(:,:,1),dx_nth(:,:,2)]'];
    end
    
    % Set callbacks for capturing next demonstration
    set(gcf,'WindowButtonDownFcn',@(h,e)button_clicked(h,e));
    set(gcf,'WindowButtonUpFcn',[]);
    set(gcf,'WindowButtonMotionFcn',[]);
    set(gcf,'Pointer','circle');
end

function ret = record_current_point(~,~)
    x = get(gca,'Currentpoint');
    x = x(1,1:2)';
    x = [x;toc];
    X = [X, x];
    hp = [hp, plot(x(1),x(2),'r.','markersize',20)];
end

end
