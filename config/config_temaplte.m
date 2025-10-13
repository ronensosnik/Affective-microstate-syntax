
# config/config_template.m (copy-paste)

```matlab
% config/config.m — user-specific paths and versions

% Project root (auto if this file is on your MATLAB path)
project_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(project_root, '..', 'code', 'matlab')));

% Data
data_root    = fullfile(project_root, '..', 'data');
data_raw     = fullfile(data_root, 'raw');          % optional
data_deriv   = fullfile(data_root, 'derivatives');  % preprocessed/epoched EEG
results_dir  = fullfile(project_root, '..', 'results');
fig_dir      = fullfile(results_dir, 'figures');
tab_dir      = fullfile(results_dir, 'tables');
stat_dir     = fullfile(results_dir, 'stats');

% Create folders if absent
cellfun(@(p) ~exist(p,'dir') && mkdir(p), {results_dir, fig_dir, tab_dir, stat_dir});

% EEGLAB/MicrostateLab
eeglab_root      = 'C:\toolboxes\eeglab2024';       % <-- EDIT
microstatelab_pl = fullfile(eeglab_root, 'plugins', 'MICROSTATELAB2.1');

% Add to path
addpath(genpath(eeglab_root));
addpath(genpath(microstatelab_pl));

% Versions (recorded in config/versions.md)
VERS.MATLAB = version();
VERS.EEGLAB = '2024.x';
VERS.MICROSTATE = '2.1';
save(fullfile(project_root, 'versions.mat'), 'VERS');
