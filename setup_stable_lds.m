% get path for setup_stable_lds
setup_stable_lds_path = which('setup_stable_lds');
setup_stable_lds_path = fileparts(setup_stable_lds_path);

% remove any auxiliary folder from the search path
restoredefaultpath();

% remove the default user-specific path
userpath('clear');

% add only the setup_stable_lds path
addpath(genpath(setup_stable_lds_path));

% 3rdParty optimizers 
% Gurobi
addpath('~/Dropbox/work/3rdParty/gurobi701/linux64/matlab/')

% Mosek
addpath('~/Dropbox/work/3rdParty/mosek/8/toolbox/');
javaaddpath('/home/medina/Dropbox/work/3rdParty/mosek/8/tools/platform/linux64x86/bin/mosekmatlab.jar','-end');