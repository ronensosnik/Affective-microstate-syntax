function run_all

% Orchestrates full analysis & figure generation

% 0) Config
run(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'config', 'config.m'));

% 1) Task performance & psychophysics
fprintf('\n[1/4] Data_analysis.m ...\n');
Data_analysis; % assumes script form; otherwise call function

% 2) EEG microstates: metrics & transitions
fprintf('\n[2/4] Microstate_analysis.m ...\n');
Microstate_analysis; % builds OUT_* and saves intermediate .mat

% 3) Tables for manuscript
fprintf('\n[3/4] Building report tables ...\n');
% Example: OUT_Coverage from workspace or .mat
% load(fullfile(stat_dir, 'OUT_Coverage.mat'), 'OUT_Coverage');
% REP = make_report_tables(OUT_Coverage);
% writetable(REP.Omnibus, fullfile(tab_dir, 'Omnibus_Coverage.csv'));
% ...

% 4) Figures (metrics + transitions)
fprintf('\n[4/4] Plotting significant metrics & transitions ...\n');
% plot_emm_significant(OUT_Occurrence, 'Occurrence', 0.05, ALL_Mean);
% GroupContrastSlides(CompStats, 'N200', {'Negative','Neutral','Positive'}, {'BP_I_Euthymic','BP_I_Depressed','BP_II_Euthymic','BP_II_Depressed','Siblings','HC'}, ...
%   'RefGroup','HC','Alpha',0.05, 'MicrostateLabels', {'A','B','C','D','E','F','G'}, ...);

fprintf('\nDone.\n');
end
