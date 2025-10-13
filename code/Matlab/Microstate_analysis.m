%% Microstate_analysis.m

% Pipeline to back-fit 7-class Custo2017 templates, segment N200/P300/LPP windows,
% and compute duration, coverage, occurrence, and directed transition counts.
% 
% Inputs (via config/config.m):
%   data_deriv : folder with preprocessed/epoched .set files organized by Group/Subject/Condition
%   results_dir: root output
% Dependencies:
%   EEGLAB, MicrostateLab (Custo2017), Statistics & ML Toolbox.
% Outputs:
%   results/stats/outputStats.mat      % per-subject microstate sequences & GFP
%   results/stats/OUT_*.mat            % OUT_Coverage / OUT_Occurrence / OUT_Duration
%   results/figures/*.png              % metrics and transition plots
% Notes:
%   - Uses continuous back-fitting (not GFP-peaks only); smoothing (lambda, b) adjustable.
%   - Windows: N200 180–300 ms; P300 300–500 ms; LPP 500–1000 ms.
%   - Replace any absolute paths with variables from config.m.

clear all
close all hidden
clc;

inputDir = 'D:\Papers\2025\In preparation\XXX (task microstates)\Data\';
saveDir = 'D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\';

ERPs           = {'N200', 'P300', 'LPP'};
Microstates  = {'Microstate_A', 'Microstate_B', 'Microstate_C', 'Microstate_D', 'Microstate_E', 'Microstate_F', 'Microstate_G'};
Groups         = {'BP_I_Depressed', 'BP_I_Euthymic', 'BP_II_Depressed', 'BP_II_Euthymic', 'HC', 'Siblings'};
Conds          = {'Negative', 'Neutral', 'Positive'};
Occurrence  = {'Occurrence_A', 'Occurrence_B', 'Occurrence_C' ,'Occurrence_D' ,'Occurrence_E' ,'Occurrence_F', 'Occurrence_G'};

N200_electrodes_names = {'Fz', 'Cz', 'C1', 'C2', 'C3', 'C4', 'FC1', 'FC2', 'FC3', 'FC4'};
P300_electrodes_names = {'AFz', 'AF3', 'AF4', 'F3', 'Fz', 'F4'};
LPP_electrodes_names =  {'CP1', 'CP2', 'CP3', 'CP4', 'P1', 'P2', 'P3', 'P4', 'PO3', 'PO4'};

ERP.window{1} = [181: 300];
ERP.window{2} = [301: 500];
ERP.window{3} = [501: 1000];

%%

%%%%%%%%%%%%%% Hierarchical clustering / averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                                                                                                              %
% Here we want to investigate topographical similarity/differences of the microstate maps themselves across subjects and groups.  %
%                                                                                                                                                                                                              %
% 1) Subject-level clustering: For each subject, pool all their trials, cluster to get subject-specific templates.                                        %
% 2) Group-level templates: Compute mean (clustering again) of subject templates within a group to get group-level templates.        %
% 3) Quantify similarity: Use similarity measures (spatial correlation) to compare subject-to-group and group-to-group templates.      %
%                                                                                                                                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ERP_to_use_for_clustering.name = '0ms_to_1200ms';
ERP_to_use_for_clustering.window = [251 1450]; % 0 to 1200ms post stimulus.

Identify_individ_microstate_maps = 0;
Identify_mean_microstate_maps_for_each_group = 0;
Identify_grand_mean_microstate_maps = 0;
Sort_grand_mean_maps = 0;
Plot_grand_mean_maps = 0;
Compare_classes_across_groups = 0;
Compare_grand_mean_Custo = 0;
Backfit_and_quantify_temp_dynam = 0;

%% Setting all parameters

groupFolders = dir(inputDir);
groupFolders = groupFolders([groupFolders.isdir]);
groupFolders = groupFolders(~matches({groupFolders.name}, {'.', '..'}));
nGroups = length(groupFolders);
groupNames = cell(1, nGroups);
dataDirs = cell(1, nGroups);

% Set clustering parameters

ClustPar.UseAAHC = false;     % true = AAHC, false = kmeans
ClustPar.MinClasses = 7;         % minimum number of clusters to identify
ClustPar.MaxClasses = 7;        % maximum number of clusters to identify
ClustPar.MaxMaps = inf;          % maximum number of data samples to use to identify clusters
ClustPar.GFPPeaks = false;    % whether clustering should be limited to global field power peaks
ClustPar.IgnorePolarity = true; % whether maps of inverted polarities should be considered part of the same cluster
ClustPar.Normalize = true;       % Set to false if using AAHC
ClustPar.Allow_early_stop = true;

ClustPar.Restarts = 100;

% Set backfitting parameters

FitPar.Classes = 7;    % cluster solutions to use for backfitting
FitPar.PeakFit = 0;    % whether to backfit only on global field power peaks
FitPar.lambda = 0.3; % smoothness penalty - ignored if FitPar.PeakFit = 1
FitPar.b = 20;            % smoothing window (ms) - ignored if FitPar.PeakFit = 1

for i = 1: nGroups
    groupDir = fullfile(inputDir, groupFolders(i).name);
    groupNames{i} = groupFolders(i).name;
    subjFolders = dir(groupDir);
    subjFolders = subjFolders(~matches({subjFolders.name}, {'.', ' ..'}));
    subjNames{i} = {subjFolders.name};
    dataDirs{i} = cellfun(@(x) fullfile(groupDir, x), subjNames{i}, 'UniformOutput', false);
end

% Start EEGLAB and find Microstates plugin files

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui');
pluginpath = fileparts(which('eegplugin_microstatelab.m'));
addpath(genpath(pluginpath));

% Load template sets

templatepath = fullfile(pluginpath, 'Templates');
Templates = dir(fullfile(templatepath, '*.set'));
MSTemplate = [];

for t = 1: numel(Templates)
    MSTemplate = eeg_store(MSTemplate, pop_loadset('filename', Templates(t).name, 'filepath', templatepath));
end

global MSTEMPLATE;
MSTEMPLATE = MSTemplate;

GroupIdx = cell(1, nGroups);
lastGroupIdx = 1;

loaded = 0;

%% Identify_individ_microstate_maps

if Identify_individ_microstate_maps

    % For each subject, load the three datasets (neutral, negative and positive valence) and update subject, group, and condition info

    for i = 1: nGroups

        for j = 1: numel(dataDirs{i})

            clear setFiles_neg setFilenames_neg setFiles_neu setFilenames_neu setFiles_pos setFilenames_pos

            setFiles_neg = dir(fullfile(dataDirs{i}{j}, '\Negative\*.set'));
            setFilenames_neg = {setFiles_neg.name};

            setFiles_neu = dir(fullfile(dataDirs{i}{j}, '\Neutral\*.set'));
            setFilenames_neu = {setFiles_neu.name};

            setFiles_pos = dir(fullfile(dataDirs{i}{j}, '\Positive\*.set'));
            setFilenames_pos = {setFiles_pos.name};

            % Load datasets

            clear EEG_neg EEG_neu EEG_pos EEG_concat

            EEG_neg = pop_loadset('filename', setFilenames_neg, 'filepath', [dataDirs{i}{j} '\Negative']);
            EEG_neg.pnts = ERP_to_use_for_clustering.window(2) - ERP_to_use_for_clustering.window(1) + 1;
            EEG_neg.xmin = ERP_to_use_for_clustering.window(1) / EEG_neg.srate;
            EEG_neg.xmax = (ERP_to_use_for_clustering.window(2) - 1) / EEG_neg.srate;
            EEG_neg.times = ERP_to_use_for_clustering.window(1): ERP_to_use_for_clustering.window(2)
            EEG_neg.data = EEG_neg.data(:, ERP_to_use_for_clustering.window(1): ERP_to_use_for_clustering.window(2), :);

            EEG_neu = pop_loadset('filename', setFilenames_neu, 'filepath', [dataDirs{i}{j} '\Neutral']);
            EEG_neu.pnts = ERP_to_use_for_clustering.window(2) - ERP_to_use_for_clustering.window(1) + 1;
            EEG_neu.xmin = ERP_to_use_for_clustering.window(1) / EEG_neu.srate;
            EEG_neu.xmax = (ERP_to_use_for_clustering.window(2) - 1) / EEG_neu.srate;
            EEG_neu.times = ERP_to_use_for_clustering.window(1): ERP_to_use_for_clustering.window(2)
            EEG_neu.data = EEG_neu.data(:, ERP_to_use_for_clustering.window(1): ERP_to_use_for_clustering.window(2), :);

            EEG_pos = pop_loadset('filename', setFilenames_pos, 'filepath', [dataDirs{i}{j} '\Positive']);
            EEG_pos.pnts = ERP_to_use_for_clustering.window(2) - ERP_to_use_for_clustering.window(1) + 1;
            EEG_pos.xmin = ERP_to_use_for_clustering.window(1) / EEG_pos.srate;
            EEG_pos.xmax = (ERP_to_use_for_clustering.window(2) - 1) / EEG_pos.srate;
            EEG_pos.times = ERP_to_use_for_clustering.window(1): ERP_to_use_for_clustering.window(2)
            EEG_pos.data = EEG_pos.data(:, ERP_to_use_for_clustering.window(1): ERP_to_use_for_clustering.window(2), :);

            EEG_concat = pop_mergeset(EEG_neg, EEG_neu);
            EEG_concat = pop_mergeset(EEG_concat, EEG_pos);
            EEG_concat.condition = 'Merged_negative_neutral_negative_conditions';
            EEG_concat.number_trials_neg_neu_pos = [EEG_neg.trials EEG_neu.trials EEG_pos.trials];

            [ALLEEG, EEG_concat, CURRENTSET] = pop_newset(ALLEEG, EEG_concat, CURRENTSET);
            currGroupIdx = lastGroupIdx: numel(ALLEEG);

            clear ppp; ppp = findstr(ALLEEG(CURRENTSET).filename, '_');

            ALLEEG(CURRENTSET).filename(ppp(end) + 1: end) = ' ';
            ALLEEG(CURRENTSET).filename(ppp(end) + 1: ppp(end) + 14) = 'all_conditions'
            ALLEEG(CURRENTSET).setname = 'Merged_negative_neutral_negative_conditions';

            % Update group and condition info for all sets

            fprintf('Updating group and condition information for group %s, subject %s...\n', groupNames{i}, subjNames{i}{j});

            for k = 1: numel(currGroupIdx)

                [EEG_concat, ALLEEG, CURRENTSET] = eeg_retrieve(ALLEEG, currGroupIdx(k));
                filename = EEG_concat.filename(1: strfind(EEG_concat.filename, '.') - 1);
                idx = strfind(filename, '_');

                if isempty(idx)
                    EEG_concat.subject = filename;
                else
                    EEG_concat.subject = filename(1: idx(2) - 1);
                end

                EEG_concat.group = groupNames{i};
                [ALLEEG, EEG_concat, CURRENTSET] = eeg_store(ALLEEG, EEG_concat, CURRENTSET);
            end

            GroupIdx{i} = [GroupIdx{i} currGroupIdx];
            lastGroupIdx = numel(ALLEEG) + 1;
        end
    end

    AllSubjects = 1: numel(ALLEEG);

    % Identify individual microstate maps for each subject

    disp('Identifying microstates for all sets...');

    [EEG, CURRENTSET] = pop_FindMSMaps(ALLEEG, AllSubjects, 'ClustPar', ClustPar);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

    Number_of_subjects = CURRENTSET(end);

    if ClustPar.MinClasses == 3
        save([saveDir ['3_15_clusters_' ERP_to_use_for_clustering.name '.mat']], 'ALLEEG', 'EEG', 'CURRENTSET', 'Number_of_subjects', 'GroupIdx', '-v7.3');
    elseif ClustPar.MinClasses == 7
        save([saveDir ['7_clusters_' ERP_to_use_for_clustering.name '.mat']], 'ALLEEG', 'EEG', 'CURRENTSET', 'Number_of_subjects', 'GroupIdx', '-v7.3');
    end

    loaded = 1;

    if ClustPar.MinClasses == 3

        threshold = 0.01; % 1% increment

        figure(1); set(gcf, 'Color', [1 1 1]);

        Index_BP_I_Depressed = 0; Index_BP_I_Euthymic = 0;
        Index_BP_II_Depressed = 0; Index_BP_II_Euthymic = 0;
        Index_Siblings = 0; Index_HC = 0;

        for subject = 1: Number_of_subjects
            if ~isempty(findstr(ALLEEG(subject).filepath, 'BP_I_Depressed'))

                Index_BP_I_Depressed = Index_BP_I_Depressed + 1;
                GEV.BP_I_Depressed(Index_BP_I_Depressed, 1) = 0/0;
                GEV.BP_I_Depressed(Index_BP_I_Depressed, 2) = 0/0;

                for pop = ClustPar.MinClasses: ClustPar.MaxClasses
                    GEV.BP_I_Depressed(Index_BP_I_Depressed, pop) = sum(ALLEEG(subject).msinfo.MSMaps(pop).ExpVar);
                end

            elseif ~isempty(findstr(ALLEEG(subject).filepath, 'BP_II_Depressed'))

                Index_BP_II_Depressed = Index_BP_II_Depressed + 1;
                GEV.BP_II_Depressed(Index_BP_II_Depressed, 1) = 0/0;
                GEV.BP_II_Depressed(Index_BP_II_Depressed, 2) = 0/0;

                for pop = ClustPar.MinClasses: ClustPar.MaxClasses
                    GEV.BP_II_Depressed(Index_BP_II_Depressed, pop) = sum(ALLEEG(subject).msinfo.MSMaps(pop).ExpVar);
                end

            elseif ~isempty(findstr(ALLEEG(subject).filepath, 'BP_I_Euthymic'))

                Index_BP_I_Euthymic = Index_BP_I_Euthymic + 1;
                GEV.BP_I_Euthymic(Index_BP_I_Euthymic, 1) = 0/0;
                GEV.BP_I_Euthymic(Index_BP_I_Euthymic, 2) = 0/0;

                for pop = ClustPar.MinClasses: ClustPar.MaxClasses
                    GEV.BP_I_Euthymic(Index_BP_I_Euthymic, pop) = sum(ALLEEG(subject).msinfo.MSMaps(pop).ExpVar);
                end

            elseif ~isempty(findstr(ALLEEG(subject).filepath, 'BP_II_Euthymic'))

                Index_BP_II_Euthymic = Index_BP_II_Euthymic + 1;
                GEV.BP_II_Euthymic(Index_BP_II_Euthymic, 1) = 0/0;
                GEV.BP_II_Euthymic(Index_BP_II_Euthymic, 2) = 0/0;

                for pop = ClustPar.MinClasses: ClustPar.MaxClasses
                    GEV.BP_II_Euthymic(Index_BP_II_Euthymic, pop) = sum(ALLEEG(subject).msinfo.MSMaps(pop).ExpVar);
                end

            elseif ~isempty(findstr(ALLEEG(subject).filepath, 'Siblings'))

                Index_Siblings = Index_Siblings + 1;
                GEV.Siblings(Index_Siblings, 1) = 0/0;
                GEV.Siblings(Index_Siblings, 2) = 0/0;

                for pop = ClustPar.MinClasses: ClustPar.MaxClasses
                    GEV.Siblings(Index_Siblings, pop) = sum(ALLEEG(subject).msinfo.MSMaps(pop).ExpVar);
                end

            elseif ~isempty(findstr(ALLEEG(subject).filepath, 'HC'))

                Index_HC = Index_HC + 1;
                GEV.HC(Index_HC, 1) = 0/0;
                GEV.HC(Index_HC, 2) = 0/0;

                for pop = ClustPar.MinClasses: ClustPar.MaxClasses
                    GEV.HC(Index_HC, pop) = sum(ALLEEG(subject).msinfo.MSMaps(pop).ExpVar);
                end

            end

        end

        clear Mean_BP_I_Depressed

        Mean_BP_I_Depressed = mean(GEV.BP_I_Depressed, 1);
        std_BP_I_Depressed = std(GEV.BP_I_Depressed, 1);

        for n = ClustPar.MinClasses + 1: ClustPar.MaxClasses
            disp(['BP_I_Depressed. ' num2str(n - 1) ' vs. ' num2str(n) ': ' num2str(round(((Mean_BP_I_Depressed(n) - Mean_BP_I_Depressed(n - 1)) / Mean_BP_I_Depressed(n - 1)) * 1000) / 10) '%']);
        end

        optimalTemplates_BP_I_Depressed = find(diff(mean(GEV.BP_I_Depressed, 1)) < threshold, 1, 'first') + 1;
        disp(['OptimalTemplates_BP_I_Depressed = ' num2str(optimalTemplates_BP_I_Depressed)]);

        disp(' ');
        disp(' ');

        clear Mean_BP_I_Euthymic

        Mean_BP_I_Euthymic = mean(GEV.BP_I_Euthymic, 1);
        std_BP_I_Euthymic = std(GEV.BP_I_Euthymic, 1);

        for n = ClustPar.MinClasses + 1: ClustPar.MaxClasses
            disp(['BP_I_Euthymic. ' num2str(n - 1) ' vs. ' num2str(n) ': ' num2str(round(((Mean_BP_I_Euthymic(n) - Mean_BP_I_Euthymic(n - 1)) / Mean_BP_I_Euthymic(n - 1)) * 1000) / 10) '%']);
        end

        optimalTemplates_BP_I_Euthymic = find(diff(mean(GEV.BP_I_Euthymic, 1)) < threshold, 1, 'first') + 1;
        disp(['optimalTemplates_BP_I_Euthymic = ' num2str(optimalTemplates_BP_I_Euthymic)]);

        disp(' ');
        disp(' ');

        clear Mean_BP_II_Depressed

        Mean_BP_II_Depressed = mean(GEV.BP_II_Depressed, 1);
        std_BP_II_Depressed = std(GEV.BP_II_Depressed, 1);

        for n = ClustPar.MinClasses + 1: ClustPar.MaxClasses
            disp(['BP_II_Depressed. ' num2str(n - 1) ' vs. ' num2str(n) ': ' num2str(round(((Mean_BP_II_Depressed(n) - Mean_BP_II_Depressed(n - 1)) / Mean_BP_II_Depressed(n - 1)) * 1000) / 10) '%']);
        end

        optimalTemplates_BP_II_Depressed = find(diff(mean(GEV.BP_II_Depressed, 1)) < threshold, 1, 'first') + 1;
        disp(['optimalTemplates_BP_II_Depressed = ' num2str(optimalTemplates_BP_II_Depressed)]);

        disp(' ');
        disp(' ');

        clear Mean_BP_II_Euthymic

        Mean_BP_II_Euthymic = mean(GEV.BP_II_Euthymic, 1);
        std_BP_II_Euthymic = std(GEV.BP_II_Euthymic, 1);

        for n = ClustPar.MinClasses + 1: ClustPar.MaxClasses
            disp(['BP_II_Euthymic. ' num2str(n - 1) ' vs. ' num2str(n) ': ' num2str(round(((Mean_BP_II_Euthymic(n) - Mean_BP_II_Euthymic(n - 1)) / Mean_BP_II_Euthymic(n - 1)) * 1000) / 10) '%']);
        end

        optimalTemplates_BP_II_Euthymic = find(diff(mean(GEV.BP_II_Euthymic)) < threshold, 1, 'first') + 1;
        disp(['optimalTemplates_BP_II_Euthymic = ' num2str(optimalTemplates_BP_II_Euthymic)]);

        disp(' ');
        disp(' ');

        clear Mean_HC
        Mean_HC = mean(GEV.HC, 1);
        std_HC = std(GEV.HC, 1);

        for n = ClustPar.MinClasses + 1: ClustPar.MaxClasses
            disp(['HC. ' num2str(n - 1) ' vs. ' num2str(n) ': ' num2str(round(((Mean_HC(n) - Mean_HC(n - 1)) / Mean_HC(n - 1)) * 1000) / 10) '%']);
        end

        optimalTemplates_HC = find(diff(mean(GEV.HC)) < threshold, 1, 'first') + 1;
        disp(['optimalTemplates_HC = ' num2str(optimalTemplates_HC)]);

        disp(' ');
        disp(' ');

        clear Mean_Siblings
        Mean_Siblings = mean(GEV.Siblings, 1);
        std_Siblings = std(GEV.Siblings, 1);

        for n = ClustPar.MinClasses + 1: ClustPar.MaxClasses
            disp(['Siblings. ' num2str(n - 1) ' vs. ' num2str(n) ': ' num2str(round(((Mean_Siblings(n) - Mean_Siblings(n - 1)) / Mean_Siblings(n - 1)) * 1000) / 10) '%']);
        end

        optimalTemplates_Siblings = find(diff(mean(GEV.Siblings)) < threshold, 1, 'first') + 1;
        disp(['optimalTemplates_Siblings = ' num2str(optimalTemplates_Siblings)]);

        disp(' ');
        disp(' ');
        disp(' ');

        disp(['7 microstate templates: BP_I_Depressed GEV = ' num2str(round(mean(GEV.BP_I_Depressed(:, 7)) * 100)) '% ' char(177) ' ' num2str(round(std(GEV.BP_I_Depressed(:, 7)) * 100)) '%']);
        disp(['7 microstate templates: BP_I_Euthymic GEV = ' num2str(round(mean(GEV.BP_I_Euthymic(:, 7)) * 100)) '% ' char(177) ' ' num2str(round(std(GEV.BP_I_Euthymic(:, 7)) * 100)) '%']);
        disp(['7 microstate templates: BP_II_Depressed GEV = ' num2str(round(mean(GEV.BP_II_Depressed(:, 7)) * 100)) '% ' char(177) ' ' num2str(round(std(GEV.BP_II_Depressed(:, 7)) * 100)) '%']);
        disp(['7 microstate templates: BP_II_Euthymic GEV = ' num2str(round(mean(GEV.BP_II_Euthymic(:, 7)) * 100)) '% ' char(177) ' ' num2str(round(std(GEV.BP_II_Euthymic(:, 7)) * 100)) '%']);
        disp(['7 microstate templates: HC GEV = ' num2str(round(mean(GEV.HC(:, 7)) * 100)) '% ' char(177) ' ' num2str(round(std(GEV.HC(:, 7)) * 100)) '%']);
        disp(['7 microstate templates: Siblings GEV = ' num2str(round(mean(GEV.Siblings(:, 7)) * 100)) '% ' char(177) ' ' num2str(round(std(GEV.Siblings(:, 7)) * 100)) '%']);

        max_GEV = ceil(max([Mean_BP_I_Depressed + std_BP_I_Depressed Mean_BP_II_Depressed + std_BP_II_Depressed Mean_BP_I_Euthymic + std_BP_I_Euthymic Mean_BP_II_Euthymic + std_BP_II_Euthymic  Mean_HC + std_HC Mean_Siblings + std_Siblings]) * 10) * 10;
        min_GEV = floor(min([Mean_BP_I_Depressed - std_BP_I_Depressed Mean_BP_II_Depressed - std_BP_II_Depressed Mean_BP_I_Euthymic - std_BP_I_Euthymic Mean_BP_II_Euthymic - std_BP_II_Euthymic  Mean_HC - std_HC Mean_Siblings - std_Siblings]) * 10) * 10;

        subplot(3, 2, 1); shadedErrorBar(1: length(Mean_BP_I_Depressed), Mean_BP_I_Depressed * 100, std_BP_I_Depressed * 100); set(gca,'XLim', [ClustPar.MinClasses ClustPar.MaxClasses], 'YLim', [min_GEV max_GEV], 'box', 'off'); title('BP\_I\_Depressed', 'FontName', 'Times_New_Roman', 'FontSize', 14);
        subplot(3, 2, 2); shadedErrorBar(1: length(Mean_BP_I_Euthymic), Mean_BP_I_Euthymic * 100, std_BP_I_Euthymic * 100); set(gca,'XLim', [ClustPar.MinClasses ClustPar.MaxClasses], 'YLim', [min_GEV max_GEV], 'box', 'off'); title('BP\_I\_Euthymic', 'FontName', 'Times_New_Roman', 'FontSize', 14);
        subplot(3, 2, 3); shadedErrorBar(1: length(Mean_BP_II_Depressed), Mean_BP_II_Depressed * 100, std_BP_II_Depressed * 100); set(gca,'XLim', [ClustPar.MinClasses ClustPar.MaxClasses], 'YLim', [min_GEV max_GEV], 'box', 'off'); title('BP\_II\_Depressed', 'FontName', 'Times_New_Roman', 'FontSize', 14);
        subplot(3, 2, 4); shadedErrorBar(1: length(Mean_BP_II_Euthymic), Mean_BP_II_Euthymic * 100, std_BP_II_Euthymic * 100); set(gca,'XLim', [ClustPar.MinClasses ClustPar.MaxClasses], 'YLim', [min_GEV max_GEV], 'box', 'off'); title('BP\_II\_Euthymic', 'FontName', 'Times_New_Roman', 'FontSize', 14);
        subplot(3, 2, 5); shadedErrorBar(1: length(Mean_HC), Mean_HC * 100, std_HC * 100); set(gca,'XLim', [ClustPar.MinClasses ClustPar.MaxClasses], 'YLim', [min_GEV max_GEV], 'box', 'off'); title('HC', 'FontName', 'Times_New_Roman', 'FontSize', 14);
        subplot(3, 2, 6); shadedErrorBar(1: length(Mean_Siblings), Mean_Siblings * 100, std_Siblings * 100); set(gca,'XLim', [ClustPar.MinClasses ClustPar.MaxClasses], 'YLim', [min_GEV max_GEV], 'box', 'off'); title('Siblings', 'FontName', 'Times_New_Roman', 'FontSize', 14);

        set(gcf, 'Color', [1 1 1]);
    end
end

%% Identify_mean_microstate_maps_for_each_group

if Identify_mean_microstate_maps_for_each_group

    if ~loaded
        load([saveDir ['\7_clusters_' ERP_to_use_for_clustering.name '.mat']]);
        loaded = 1;
    end

    GroupMeanIdx = [];
    CURRENTSET = Number_of_subjects;

    for i = 1: nGroups

        fprintf('Identifying mean maps for group %s\n', groupNames{i});

        EEG = pop_CombMSMaps(ALLEEG, GroupIdx{i}, 'MeanName', ['Mean_' groupNames{i}], 'IgnorePolarity', ClustPar.IgnorePolarity);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
        GroupMeanIdx = [GroupMeanIdx CURRENTSET];
    end

    save([saveDir ['\7_clusters_' ERP_to_use_for_clustering.name '.mat']], 'ALLEEG', 'EEG', 'CURRENTSET', 'GroupIdx', 'GroupMeanIdx', 'Number_of_subjects', '-v7.3');
end

%% Identify_grand_mean_microstate_maps

if Identify_grand_mean_microstate_maps

    if ~loaded
        load([saveDir ['\7_clusters_' ERP_to_use_for_clustering.name '.mat']]);
        loaded = 1;
    end

    CURRENTSET = Number_of_subjects + 6;

    fprintf('Identifying grand mean maps');

    EEG = pop_CombMSMaps(ALLEEG, GroupMeanIdx, 'MeanName', 'Grand_mean', 'IgnorePolarity', ClustPar.IgnorePolarity);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');

    save([saveDir ['\7_clusters_' ERP_to_use_for_clustering.name '.mat']], 'ALLEEG', 'EEG', 'CURRENTSET', 'GroupIdx', 'GroupMeanIdx', 'Number_of_subjects', '-v7.3');
end

%% Sort grand mean maps by the 2017 Custo maps and then sort all individual maps by the grand mean map

if Sort_grand_mean_maps

    if ~loaded
        load([saveDir ['\7_clusters_' ERP_to_use_for_clustering.name '.mat']]);
        loaded = 1;
    end

    % Sort the mean maps by the specified published template(s)
    % First, sort the 7-class solution of the grand mean by the 2017 Custo maps.

    disp(['Grand mean (across groups) mean (across seven classes) shared variance: ' num2str(round(mean(ALLEEG(CURRENTSET).msinfo.MSMaps(7).SharedVar) * 10000) / 100)]);

    [ALLEEG, EEG, CURRENTSET] = pop_SortMSMaps(ALLEEG, CURRENTSET, 'TemplateSet', 'Custo2017', 'Classes', 7, 'IgnorePolarity', 1);

    % Next, sort the seven classes of each individual template in each group by the seven classes of the grand mean

    [ALLEEG, EEG, CURRENTSET] = pop_SortMSMaps(ALLEEG, 1: CURRENTSET - 1, 'TemplateSet', CURRENTSET, 'Classes', 7, 'IgnorePolarity', 1);

    save ([saveDir ['\7_clusters_' ERP_to_use_for_clustering.name '.mat']], 'ALLEEG', 'EEG', 'CURRENTSET', 'GroupIdx', 'GroupMeanIdx', 'Number_of_subjects',  '-v7.3');

end

%%

if Plot_grand_mean_maps

    eeglab;

    if ~loaded
        load([saveDir ['\7_clusters_' ERP_to_use_for_clustering.name '.mat']]);
        loaded = 1;
    end

    pop_ShowIndMSMaps(ALLEEG, Number_of_subjects + 1: Number_of_subjects + 7, 'Classes', 7);

    load('D:\eeglab2025.0.0\plugins\MICROSTATELAB2.1\Templates\Custo.mat');
    pop_ShowIndMSMaps(Custo, 1, 'Classes', 7);
end

%% Compare classes across groups

if Compare_classes_across_groups
    eeglab
    if ~loaded
        load([saveDir ['\7_clusters_' ERP_to_use_for_clustering.name '.mat']]);
        loaded = 1;
    end

    % Compare group means

    sharedVarTable = pop_CompareMSMaps(ALLEEG, 'MeanSets', Number_of_subjects + 1: Number_of_subjects + 6, 'Classes', 7, 'gui', 0);
    sharedVarTable = table2array(sharedVarTable);

    disp(' ');
    disp('Shared variance - compare all group means. Congruent')
    disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
    disp(' ');

    A_microstate_mean = mean([sharedVarTable(1, 8) sharedVarTable(1, 15) sharedVarTable(1, 22)   sharedVarTable(1, 29)   sharedVarTable(1, 36) ...
        sharedVarTable(8, 15) sharedVarTable(8, 22)   sharedVarTable(8, 29)   sharedVarTable(8, 36) ...
        sharedVarTable(15, 22)  sharedVarTable(15, 29) sharedVarTable(15, 36) ...
        sharedVarTable(22, 29) sharedVarTable(22, 36) ...
        sharedVarTable(29, 36)]);

    A_microstate_std = std([sharedVarTable(1, 8) sharedVarTable(1, 15) sharedVarTable(1, 22)   sharedVarTable(1, 29)   sharedVarTable(1, 36) ...
        sharedVarTable(8, 15) sharedVarTable(8, 22)   sharedVarTable(8, 29)   sharedVarTable(8, 36) ...
        sharedVarTable(15, 22)  sharedVarTable(15, 29) sharedVarTable(15, 36) ...
        sharedVarTable(22, 29) sharedVarTable(22, 36) ...
        sharedVarTable(29, 36)]);

    B_microstate_mean = mean([sharedVarTable(2, 9) sharedVarTable(2, 16) sharedVarTable(2, 23)   sharedVarTable(2, 30)   sharedVarTable(2, 37) ...
        sharedVarTable(9, 16) sharedVarTable(9, 23)   sharedVarTable(9, 30)   sharedVarTable(9, 37) ...
        sharedVarTable(16, 23)  sharedVarTable(16, 30) sharedVarTable(16, 37) ...
        sharedVarTable(23, 30) sharedVarTable(23, 37) ...
        sharedVarTable(30, 37)]);

    B_microstate_std = std([sharedVarTable(2, 9) sharedVarTable(2, 16) sharedVarTable(2, 23)   sharedVarTable(2, 30)   sharedVarTable(2, 37) ...
        sharedVarTable(9, 16) sharedVarTable(9, 23)   sharedVarTable(9, 30)   sharedVarTable(9, 37) ...
        sharedVarTable(16, 23)  sharedVarTable(16, 30) sharedVarTable(16, 37) ...
        sharedVarTable(23, 30) sharedVarTable(23, 37) ...
        sharedVarTable(30, 37)]);

    C_microstate_mean = mean([sharedVarTable(3, 10) sharedVarTable(3, 17) sharedVarTable(3, 24)    sharedVarTable(3, 31)    sharedVarTable(3, 38) ...
        sharedVarTable(10, 17) sharedVarTable(10, 24)  sharedVarTable(10, 31)  sharedVarTable(10, 38) ...
        sharedVarTable(17, 24)  sharedVarTable(17, 31)   sharedVarTable(17, 38) ...
        sharedVarTable(24, 31)   sharedVarTable(24, 38) ...
        sharedVarTable(31, 38)]);

    C_microstate_std = std([sharedVarTable(3, 10) sharedVarTable(3, 17) sharedVarTable(3, 24)    sharedVarTable(3, 31)    sharedVarTable(3, 38) ...
        sharedVarTable(10, 17) sharedVarTable(10, 24)  sharedVarTable(10, 31)  sharedVarTable(10, 38) ...
        sharedVarTable(17, 24)  sharedVarTable(17, 31)   sharedVarTable(17, 38) ...
        sharedVarTable(24, 31)   sharedVarTable(24, 38) ...
        sharedVarTable(31, 38)]);

    D_microstate_mean = mean([sharedVarTable(4, 11) sharedVarTable(4, 18) sharedVarTable(4, 25)    sharedVarTable(4, 32)    sharedVarTable(4, 39) ...
        sharedVarTable(11, 18) sharedVarTable(11, 25)  sharedVarTable(11, 32)  sharedVarTable(11, 39) ...
        sharedVarTable(18, 25)  sharedVarTable(18, 32)   sharedVarTable(18, 39) ...
        sharedVarTable(25, 32)   sharedVarTable(25, 39) ...
        sharedVarTable(32, 39)]);

    D_microstate_std = std([sharedVarTable(4, 11) sharedVarTable(4, 18) sharedVarTable(4, 25)    sharedVarTable(4, 32)    sharedVarTable(4, 39) ...
        sharedVarTable(11, 18) sharedVarTable(11, 25)  sharedVarTable(11, 32)  sharedVarTable(11, 39) ...
        sharedVarTable(18, 25)  sharedVarTable(18, 32)   sharedVarTable(18, 39) ...
        sharedVarTable(25, 32)   sharedVarTable(25, 39) ...
        sharedVarTable(32, 39)]);

    E_microstate_mean = mean([sharedVarTable(5, 12) sharedVarTable(5, 19) sharedVarTable(5, 26)    sharedVarTable(5, 33)    sharedVarTable(5, 40) ...
        sharedVarTable(12, 19) sharedVarTable(12, 26)  sharedVarTable(12, 33)  sharedVarTable(12, 40) ...
        sharedVarTable(19, 26)  sharedVarTable(19, 33)   sharedVarTable(19, 40) ...
        sharedVarTable(26, 33)   sharedVarTable(26, 40) ...
        sharedVarTable(33, 40)]);

    E_microstate_std = std([sharedVarTable(5, 12) sharedVarTable(5, 19) sharedVarTable(5, 26)    sharedVarTable(5, 33)    sharedVarTable(5, 40) ...
        sharedVarTable(12, 19) sharedVarTable(12, 26)  sharedVarTable(12, 33)  sharedVarTable(12, 40) ...
        sharedVarTable(19, 26)  sharedVarTable(19, 33)   sharedVarTable(19, 40) ...
        sharedVarTable(26, 33)   sharedVarTable(26, 40) ...
        sharedVarTable(33, 40)]);

    F_microstate_mean = mean([sharedVarTable(6, 13) sharedVarTable(6, 20) sharedVarTable(6, 27)    sharedVarTable(6, 34)    sharedVarTable(6, 41) ...
        sharedVarTable(13, 20) sharedVarTable(13, 27)  sharedVarTable(13, 34)  sharedVarTable(13, 41) ...
        sharedVarTable(20, 27)  sharedVarTable(20, 34)   sharedVarTable(20, 41) ...
        sharedVarTable(27, 34)   sharedVarTable(27, 41) ...
        sharedVarTable(34, 41)]);

    F_microstate_std = std([sharedVarTable(6, 13) sharedVarTable(6, 20) sharedVarTable(6, 27)    sharedVarTable(6, 34)    sharedVarTable(6, 41) ...
        sharedVarTable(13, 20) sharedVarTable(13, 27)  sharedVarTable(13, 34)  sharedVarTable(13, 41) ...
        sharedVarTable(20, 27)  sharedVarTable(20, 34)   sharedVarTable(20, 41) ...
        sharedVarTable(27, 34)   sharedVarTable(27, 41) ...
        sharedVarTable(34, 41)]);

    G_microstate_mean = mean([sharedVarTable(7, 14) sharedVarTable(7, 21) sharedVarTable(7, 28)    sharedVarTable(7, 35)    sharedVarTable(7, 42) ...
        sharedVarTable(14, 21) sharedVarTable(14, 28)  sharedVarTable(14, 35)  sharedVarTable(14, 42) ...
        sharedVarTable(21, 28)  sharedVarTable(21, 35)   sharedVarTable(21, 42) ...
        sharedVarTable(28, 35)   sharedVarTable(28, 42) ...
        sharedVarTable(35, 42)]);

    G_microstate_std = std([sharedVarTable(7, 14) sharedVarTable(7, 21) sharedVarTable(7, 28)    sharedVarTable(7, 35)    sharedVarTable(7, 42) ...
        sharedVarTable(14, 21) sharedVarTable(14, 28)  sharedVarTable(14, 35)  sharedVarTable(14, 42) ...
        sharedVarTable(21, 28)  sharedVarTable(21, 35)   sharedVarTable(21, 42) ...
        sharedVarTable(28, 35)   sharedVarTable(28, 42) ...
        sharedVarTable(35, 42)]);

    disp(['A_microstate = ' num2str(A_microstate_mean) ' ' char(177) ' ' num2str(A_microstate_std)]);
    disp(['B_microstate = ' num2str(B_microstate_mean) ' ' char(177) ' ' num2str(B_microstate_std)]);
    disp(['C_microstate = ' num2str(C_microstate_mean) ' ' char(177) ' ' num2str(C_microstate_std)]);
    disp(['D_microstate = ' num2str(D_microstate_mean) ' ' char(177) ' ' num2str(D_microstate_std)]);
    disp(['E_microstate = ' num2str(E_microstate_mean) ' ' char(177) ' ' num2str(E_microstate_std)]);
    disp(['F_microstate = ' num2str(F_microstate_mean) ' ' char(177) ' ' num2str(F_microstate_std)]);
    disp(['G_microstate = ' num2str(G_microstate_mean) ' ' char(177) ' ' num2str(G_microstate_std)]);

    disp(' ');
    disp(' ');
    disp(' ');

end

%% Compare grand mean vs. Custo

if Compare_grand_mean_Custo

    eeglab

    if ~loaded
        load([saveDir ['\7_clusters_' ERP_to_use_for_clustering.name '.mat']]);
        loaded = 1;
    end

    % Compare grand mean and Custo et al template

    sharedVarTable_GrandMean = pop_CompareMSMaps(ALLEEG, 'MeanSets', 72, 'PublishedSets', 'Custo2017', 'Classes', 7, 'Filename', [saveDir '\sharedvars.mat'], 'gui', 0);
    sharedVarTable_GrandMean = table2array(sharedVarTable_GrandMean);

    % Compute similarity index (shared variance across the seven global microstates of each group)

    Similarity_incongruent_microstate = []; Similarity_A_microstate = []; Similarity_B_microstate = []; Similarity_C_microstate = []; Similarity_D_microstate = []; Similarity_E_microstate = []; Similarity_F_microstate = []; Similarity_G_microstate = [];

    Similarity_A_microstate = sharedVarTable_GrandMean(1, 8);
    Similarity_B_microstate = sharedVarTable_GrandMean(2, 9);
    Similarity_C_microstate = sharedVarTable_GrandMean(3, 10);
    Similarity_D_microstate = sharedVarTable_GrandMean(4, 11);
    Similarity_E_microstate = sharedVarTable_GrandMean(5, 12);
    Similarity_F_microstate = sharedVarTable_GrandMean(6, 13);
    Similarity_G_microstate = sharedVarTable_GrandMean(7, 14);

    for row = 1: 13
        for col = row + 1: 7: 14
            try
                Similarity_incongruent_microstate = [Similarity_incongruent_microstate sharedVarTable_GrandMean(row, col: col + 5)];
            catch
                Similarity_incongruent_microstate = [Similarity_incongruent_microstate sharedVarTable_GrandMean(row, col: 14)];
            end
        end
    end

    disp('Shared variance - grand mean vs. Custo template. Congruent')
    disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ')
    disp(' ');

    disp(['A_microstate = ' num2str(mean(Similarity_A_microstate))]);
    disp(['B_microstate = ' num2str(mean(Similarity_B_microstate))]);
    disp(['C_microstate = ' num2str(mean(Similarity_C_microstate))]);
    disp(['D_microstate = ' num2str(mean(Similarity_D_microstate))]);
    disp(['E_microstate = ' num2str(mean(Similarity_E_microstate))]);
    disp(['F_microstate = ' num2str(mean(Similarity_F_microstate))]);
    disp(['G_microstate = ' num2str(mean(Similarity_G_microstate))]);

    disp(' ');
    disp(' ');
    disp(' ');

    disp('Shared variance - grand mean vs. Custo template. Incongruent')
    disp(' - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - ')

    disp(['All_microstates = ' num2str(mean(Similarity_incongruent_microstate))]);

    sharedVarTablee= sharedVarTable_GrandMean;

    sharedVarTable_plot(1, :) =    [sharedVarTablee(1, 1)      sharedVarTablee(1, 8)      sharedVarTablee(1, 2)       sharedVarTablee(1, 9)        sharedVarTablee(1, 3)         sharedVarTablee(1, 10)          sharedVarTablee(1, 4)          sharedVarTablee(1, 11)           sharedVarTablee(1, 5)           sharedVarTablee(1, 12)           sharedVarTablee(1, 6)             sharedVarTablee(1, 13)             sharedVarTablee(1, 7)           sharedVarTablee(1, 14)];
    sharedVarTable_plot(2, :) =    [0/0                                    sharedVarTablee(8, 8)      sharedVarTablee(8, 2)       sharedVarTablee(8, 9)         sharedVarTablee(8, 3)         sharedVarTablee(8, 10)          sharedVarTablee(8, 4)          sharedVarTablee(8, 11)           sharedVarTablee(8, 5)           sharedVarTablee(8, 12)           sharedVarTablee(8, 6)             sharedVarTablee(8, 13)             sharedVarTablee(8, 7)           sharedVarTablee(8, 14)];
    sharedVarTable_plot(3, :) =    [0/0                                    0/0                                    sharedVarTablee(2, 2)        sharedVarTablee(2, 9)        sharedVarTablee(2, 3)         sharedVarTablee(2, 10)          sharedVarTablee(2, 4)          sharedVarTablee(2, 11)           sharedVarTablee(2, 5)           sharedVarTablee(2, 12)           sharedVarTablee(2, 6)             sharedVarTablee(2, 13)             sharedVarTablee(2, 7)           sharedVarTablee(2, 14)];
    sharedVarTable_plot(4, :) =    [0/0                                    0/0                                    0/0                                      sharedVarTablee(9, 9)         sharedVarTablee(9, 3)         sharedVarTablee(9, 10)          sharedVarTablee(9, 4)          sharedVarTablee(9, 11)           sharedVarTablee(9, 5)           sharedVarTablee(9, 12)           sharedVarTablee(9, 6)             sharedVarTablee(9, 13)             sharedVarTablee(9, 7)           sharedVarTablee(9, 14)];
    sharedVarTable_plot(5, :) =    [0/0                                    0/0                                    0/0                                      0/0                                       sharedVarTablee(3, 3)         sharedVarTablee(3, 10)          sharedVarTablee(3, 4)          sharedVarTablee(3, 11)           sharedVarTablee(3, 5)           sharedVarTablee(3, 12)           sharedVarTablee(3, 6)             sharedVarTablee(3, 13)             sharedVarTablee(3, 7)           sharedVarTablee(3, 14)];
    sharedVarTable_plot(6, :) =    [0/0                                    0/0                                    0/0                                      0/0                                       0/0                                       sharedVarTablee(10, 10)        sharedVarTablee(10, 4)        sharedVarTablee(10, 11)         sharedVarTablee(10, 5)         sharedVarTablee(10, 12)         sharedVarTablee(10, 6)           sharedVarTablee(10, 13)           sharedVarTablee(10, 7)         sharedVarTablee(10, 14)];
    sharedVarTable_plot(7, :) =    [0/0                                    0/0                                    0/0                                      0/0                                       0/0                                       0/0                                           sharedVarTablee(4, 4)          sharedVarTablee(4, 11)           sharedVarTablee(4, 5)           sharedVarTablee(4, 12)           sharedVarTablee(4, 6)             sharedVarTablee(4, 13)             sharedVarTablee(4, 7)           sharedVarTablee(4, 14)];
    sharedVarTable_plot(8, :) =    [0/0                                    0/0                                    0/0                                      0/0                                       0/0                                       0/0                                           0/0                                        sharedVarTablee(11, 11)         sharedVarTablee(11, 5)         sharedVarTablee(11, 12)          sharedVarTablee(11, 6)           sharedVarTablee(11, 13)           sharedVarTablee(11, 7)         sharedVarTablee(11, 14)];
    sharedVarTable_plot(9, :) =    [0/0                                    0/0                                    0/0                                      0/0                                       0/0                                       0/0                                           0/0                                        0/0                                           sharedVarTablee(5, 5)           sharedVarTablee(5, 12)           sharedVarTablee(5, 6)             sharedVarTablee(5, 13)             sharedVarTablee(5, 7)           sharedVarTablee(5, 14)];
    sharedVarTable_plot(10, :) =  [0/0                                    0/0                                    0/0                                      0/0                                       0/0                                       0/0                                           0/0                                        0/0                                           0/0                                         sharedVarTablee(12, 12)         sharedVarTablee(12, 6)           sharedVarTablee(12, 13)           sharedVarTablee(12, 7)         sharedVarTablee(12, 14)];
    sharedVarTable_plot(11, :) =  [0/0                                    0/0                                    0/0                                      0/0                                       0/0                                       0/0                                           0/0                                        0/0                                           0/0                                          0/0                                           sharedVarTablee(6, 6)             sharedVarTablee(6, 13)             sharedVarTablee(6, 7)           sharedVarTablee(6, 14)];
    sharedVarTable_plot(12, :) =  [0/0                                    0/0                                    0/0                                      0/0                                       0/0                                       0/0                                           0/0                                        0/0                                           0/0                                          0/0                                           0/0                                           sharedVarTablee(13, 13)           sharedVarTablee(13, 7)         sharedVarTablee(13, 14)];
    sharedVarTable_plot(13, :) =  [0/0                                    0/0                                    0/0                                      0/0                                       0/0                                       0/0                                           0/0                                        0/0                                           0/0                                          0/0                                           0/0                                           0/0                                              sharedVarTablee(7, 7)           sharedVarTablee(7, 14)];
    sharedVarTable_plot(14, :) =  [0/0                                    0/0                                    0/0                                      0/0                                       0/0                                       0/0                                           0/0                                        0/0                                           0/0                                          0/0                                           0/0                                           0/0                                              0/0                                          sharedVarTablee(14, 14)];

    figure(2);

    set(gcf, 'Color', [1 1 1]);
    heatmap(sharedVarTable_plot, 'CellLabelColor', 'none', 'MissingDataColor', [1 1 1], 'GridVisible', 'off', 'ColorBarVisible' ,'off', ...
        'Title', ['Comparison. Grand mean vs. Custo template. ' ERP_to_use_for_clustering.name]);
    colormap('jet');
    colorbar;

end

if Backfit_and_quantify_temp_dynam

    if ~loaded
        load([saveDir ['\7_clusters_' ERP_to_use_for_clustering.name '.mat']]);
        loaded = 1;
    end

    Groups = {'BP_II_Depressed', 'BP_II_Euthymic', 'BP_I_Depressed', 'BP_I_Euthymic', 'HC', 'Siblings'};

    %% Backfit using Custo2017 and quantify temporal dynamics

    disp('Backfitting and extracting temporal dynamics...');

    [EEG, CURRENTSET] = pop_FitMSMaps(ALLEEG, 1: Number_of_subjects, 'TemplateSet', 'Custo2017', 'FitPar', FitPar);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

    save([saveDir ['\7_clusters_' ERP_to_use_for_clustering.name '.mat']], 'ALLEEG', 'EEG', 'CURRENTSET', 'GroupIdx', 'GroupMeanIdx', 'Number_of_subjects',  '-v7.3');

    for aaa = 1: 6
        subj_numm(aaa) = size(GroupIdx{aaa}, 2);
    end

    subj_index = cumsum(subj_numm);
    subj_index(2: end + 1) = subj_index;
    subj_index(1) = 0;

    for group = 1: 6
        inddd = 0;

        for num_subjects = subj_index(group) + 1: subj_index(group + 1)
            inddd = inddd + 1;
            clear number_trials_neg_neu_pos; number_trials_neg_neu_pos = cumsum(ALLEEG(num_subjects).number_trials_neg_neu_pos);

            outputStats.(Groups{group}).negative(inddd).MSClass = ALLEEG(num_subjects).msinfo.MSStats(7).MSClass(1: end, 1: number_trials_neg_neu_pos(1));
            outputStats.(Groups{group}).negative(inddd).GFP = ALLEEG(num_subjects).msinfo.MSStats(7).GFP(1: end, 1: number_trials_neg_neu_pos(1));

            outputStats.(Groups{group}).neutral(inddd).MSClass = ALLEEG(num_subjects).msinfo.MSStats(7).MSClass(1: end, number_trials_neg_neu_pos(1) + 1: number_trials_neg_neu_pos(2));
            outputStats.(Groups{group}).neutral(inddd).GFP = ALLEEG(num_subjects).msinfo.MSStats(7).GFP(1: end, number_trials_neg_neu_pos(1) + 1: number_trials_neg_neu_pos(2));

            outputStats.(Groups{group}).positive(inddd).MSClass = ALLEEG(num_subjects).msinfo.MSStats(7).MSClass(1: end, number_trials_neg_neu_pos(2) + 1: number_trials_neg_neu_pos(3));
            outputStats.(Groups{group}).positive(inddd).GFP = ALLEEG(num_subjects).msinfo.MSStats(7).GFP(1: end, number_trials_neg_neu_pos(2) + 1: number_trials_neg_neu_pos(3));
        end
    end

    save(fullfile(saveDir, ['outputStats.mat']), 'outputStats');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Within each window, labels were clipped to the window bounds.
% For each microstate class we computed:
% (i) coverage, the proportion of samples in the window carrying that label (segments spanning a boundary contributed by their overlap only);
% (ii) occurrence, the number of onsets whose first sample fell within the window, divided by the window length to yield a rate (Hz);
% (iii) mean duration, the average duration of segments whose onsets occurred within the window (segments truncated by the window boundary were excluded from the duration calculation).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('D:\Papers\2025\In preparation\XXX (task microstates)\Data\France_data.mat');
load(fullfile(saveDir, ['outputStats.mat']), 'outputStats');

% leaf template (the 5 fields at the bottom level)

leaf = struct('num_occurrences', [], 'duration', [], 'coverage', [], 'number', [], 'age', []);

% build the nested struct

Data = struct();

for i = 1: numel(ERPs)
    e = ERPs{i};
    Data.(e) = struct();

    for j = 1: numel(Microstates)
        m = Microstates{j};
        Data.(e).(m) = struct();

        for k = 1: numel(Groups)
            g = Groups{k};
            Data.(e).(m).(g) = struct();

            for t = 1: numel(Conds)
                c = Conds{t};
                Data.(e).(m).(g).(c) = repmat({leaf}, 1, 1);  % create the 5 fields here
                Mean.(e).(m).(g).(c) = repmat({leaf}, 1, 1); % create the 5 fields here
            end
        end
    end
end

for erp = 1: 3
    for micro = 1: 7
        for gro = 1: 6
            for condd = 1: 3

                clear vvv; vvv = {outputStats.(Groups{gro}).([lower(Conds{condd}(1)) Conds{condd}(2: end)]).MSClass};

                for subj_num = 1: size(vvv, 2)

                    for trial_num = 1: size(vvv{subj_num}, 2)

                        clear x d y run_starts run_ends z

                        x = vvv{subj_num}(ERP.window{erp}(1): ERP.window{erp}(end), trial_num)';
                        d = [true diff(x) ~= 0];
                        y = x(d);

                        run_starts = find(d);
                        run_ends = [run_starts(2: end) - 1, length(x)];

                        z = run_ends - run_starts + 1;

                        % occurrence

                        for pppp = 1: 7
                            aaa(pppp) = sum(y == pppp);
                        end

                        for pppp = 1: 7
                            if y(1) == pppp && vvv{subj_num}(ERP.window{erp}(1) - 1, trial_num) == pppp
                                aaa(pppp) = aaa(pppp) - 1;
                            end
                        end

                        try
                            Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.num_occurrences = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.num_occurrences aaa(micro) * (1000 / length(ERP.window{erp}))];
                        catch
                            Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.num_occurrences = [];
                            Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.num_occurrences = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.num_occurrences aaa(micro) * (1000 / length(ERP.window{erp}))];
                        end

                        % mean duration

                        clear x_ee; x_ee = x;

                        for pppp = 1: 7
                            if y(1) == pppp && vvv{subj_num}(ERP.window{erp}(1) - 1, trial_num) == pppp
                                x_ee(1: z(1)) = -9999;
                            end

                            if y(end) == pppp && vvv{subj_num}(ERP.window{erp}(end) + 1, trial_num) == pppp
                                x_ee(run_starts(end): run_ends(end)) = -9999;
                            end
                        end

                        if x_ee(1) ~= -9999
                            d_ee = [true diff(x_ee) ~= 0];
                        else
                            d_ee = [false diff(x_ee) ~= 0];
                        end

                        y_ee = x_ee(d_ee);

                        run_starts_ee = find(d_ee);
                        run_ends_ee = [run_starts_ee(2: end) - 1, length(x_ee)];

                        z_ee = run_ends_ee - run_starts_ee + 1;

                        fff = find(y_ee == -9999);

                        if ~isempty(fff)
                            y_ee(fff) = []; run_starts_ee(fff) = []; run_ends_ee(fff) = []; z_ee(fff) = [];
                        end

                        try
                            Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.duration = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.duration mean(z_ee(y_ee == micro))];
                        catch
                            Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.duration = [];
                            Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.duration = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.duration mean(z_ee(y_ee == micro))];
                        end

                        % coverage

                        try
                            Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.coverage = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.coverage (sum(z(y == micro))) / (ERP.window{erp}(end) - ERP.window{erp}(1) + 1)];
                        catch
                            Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.coverage = [];
                            Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.coverage = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.coverage (sum(z(y == micro))) / (ERP.window{erp}(end) - ERP.window{erp}(1) + 1)];
                        end

                        % subject number & age

                        switch Groups{gro}
                            case 'BP_I_Depressed'
                                try
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_I.Depressed{subj_num}.age];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_I.Depressed{subj_num}.number];
                                catch
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_I.Depressed{subj_num}.age];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_I.Depressed{subj_num}.number];
                                end

                            case 'BP_II_Depressed'
                                try
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_II.Depressed{subj_num}.age];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_II.Depressed{subj_num}.number];
                                catch
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_II.Depressed{subj_num}.age];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_II.Depressed{subj_num}.number];
                                end

                            case 'BP_I_Euthymic'
                                try
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_I.Euthymic{subj_num}.age];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_I.Euthymic{subj_num}.number];
                                catch
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_I.Euthymic{subj_num}.age];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_I.Euthymic{subj_num}.number];
                                end

                            case 'BP_II_Euthymic'
                                try
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_II.Euthymic{subj_num}.age];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_II.Euthymic{subj_num}.number];
                                catch
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_II.Euthymic{subj_num}.age];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_II.Euthymic{subj_num}.number];
                                end

                            case 'HC'
                                try
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.HC{subj_num}.age];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.HC{subj_num}.number];
                                catch
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.HC{subj_num}.age];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.HC{subj_num}.number];
                                end

                            case 'Siblings'
                                try
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.Siblings{subj_num}.age];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.Siblings{subj_num}.number];
                                catch
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.Siblings{subj_num}.age];
                                    Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.Siblings{subj_num}.number];
                                end

                        end
                    end

                end

                for subj_num = 1: size(vvv, 2)

                    Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.num_occurrences = mean(Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.num_occurrences);
                    Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.duration = nanmean(Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.duration);
                    Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.coverage = mean(Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.coverage) * 100;
                    Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.group = Groups{gro};
                    Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.condition = Conds{condd};

                    switch Groups{gro}

                        case 'BP_I_Depressed'
                            Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = France.BP_I.Depressed{subj_num}.number;
                            Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = France.BP_I.Depressed{subj_num}.age;
                        case 'BP_II_Depressed'
                            Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = France.BP_II.Depressed{subj_num}.number;
                            Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = France.BP_II.Depressed{subj_num}.age;
                        case 'BP_I_Euthymic'
                            Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = France.BP_I.Euthymic{subj_num}.number;
                            Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = France.BP_I.Euthymic{subj_num}.age;
                        case 'BP_II_Euthymic'
                            Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = France.BP_II.Euthymic{subj_num}.number;
                            Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = France.BP_II.Euthymic{subj_num}.age;
                        case 'HC'
                            Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = France.HC{subj_num}.number;
                            Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = France.HC{subj_num}.age;
                        case 'Siblings'
                            Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = France.Siblings{subj_num}.number;
                            Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = France.Siblings{subj_num}.age;

                    end

                end

            end
        end
    end
end

%% Populate variables comprising all groups and conditions data for each ERP and microstate

for erp = 1: 3

    eName = ERPs{erp};

    for micro = 1: 7

        mName = Microstates{micro};

        % init accumulators once per (ERP,Microstate)

        ALL_Mean.(eName).(mName).Age                      = [];
        ALL_Mean.(eName).(mName).Number                = [];
        ALL_Mean.(eName).(mName).Coverage             = [];
        ALL_Mean.(eName).(mName).Duration               = [];
        ALL_Mean.(eName).(mName).Num_occurrence = [];
        ALL_Mean.(eName).(mName).Group                  = {};
        ALL_Mean.(eName).(mName).Condition             = {};

        for gro = 1: 6
            gName = Groups{gro};

            for condd = 1: 3
                cName = Conds{condd};

                % this should be a 1xN cell array of subject structs

                S = Mean.(eName).(mName).(gName).(cName);

                if isempty(S)
                    continue;
                end

                % extract vectors from the cell array of structs

                ages = cellfun(@(s) s.age, S);
                nums = cellfun(@(s) s.number, S);
                covs = cellfun(@(s) s.coverage, S);
                durs = cellfun(@(s) s.duration, S);
                occs = cellfun(@(s) s.num_occurrences, S);
                grps = cellfun(@(s) s.group, S, 'UniformOutput', false);
                conds = cellfun(@(s) s.condition, S, 'UniformOutput', false);

                % append (concatenate) to accumulators

                ALL_Mean.(eName).(mName).Age                     = [ALL_Mean.(eName).(mName).Age,            ages];
                ALL_Mean.(eName).(mName).Number               = [ALL_Mean.(eName).(mName).Number,         nums];
                ALL_Mean.(eName).(mName).Coverage            = [ALL_Mean.(eName).(mName).Coverage,       covs];
                ALL_Mean.(eName).(mName).Duration              = [ALL_Mean.(eName).(mName).Duration,       durs];
                ALL_Mean.(eName).(mName).Num_occurrence = [ALL_Mean.(eName).(mName).Num_occurrence, occs];
                ALL_Mean.(eName).(mName).Group                  = [ALL_Mean.(eName).(mName).Group,          grps];
                ALL_Mean.(eName).(mName).Condition             = [ALL_Mean.(eName).(mName).Condition,      conds];
            end
        end
    end
end

save('D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\Metrics_ALL.mat', 'ALL_Mean');

% For each of the three measures, we run 21 sets of models (7 microstates * 3 ERPs)

Out_Duration = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Duration', 0.05);
Out_Occurrence = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Occurrence', 0.05);
Out_Coverage = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Coverage', 0.05);

% Report

print_summary_measure(Out_Duration, 0.05);
print_summary_measure(Out_Occurrence, 0.05);
print_summary_measure(Out_Coverage, 0.05);

% Generate plots for any FDR-significant effects

plot_emm_significant(Out_Duration,      'Duration',      0.05, ALL_Mean);
plot_emm_significant(Out_Occurrence, 'Occurrence', 0.05, ALL_Mean);
plot_emm_significant(Out_Coverage,    'Coverage',    0.05, ALL_Mean);

%% ======================================================================= %%
%  TASK MICROSTATES — TWO TRANSITION ANALYSES (AGE-ADJUSTED, FDR CONTROL)
%  Q1. Within each cell (ERP×Group×Condition), which X->Y transitions
%         deviate from independence?  [Poisson GLMs with independence offset]
%  Q2. Anchored contrasts: For each ERP×Condition, which X->Y transitions in
%         each target group differ from a chosen reference (age-anchored)?
%
%  Outputs:
%   - TransStats_cellFDR: per-cell deviations (q-values, edges for plots)
%   - CompStats: Target vs Reference anchored contrasts (q-values, edges)
%
%% ======================================================================== %%

%% ------------------------------ I/O ---------------------------------- %%
data_path = 'D:\Papers\2025\In preparation\XXX (task microstates)\Data\France_data.mat';
stats_path = 'D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\outputStats.mat';
save_path = 'D:\Papers\2025\In preparation\XXX (task microstates)\Analysis';

load(data_path, 'France'); % ages & metadata (struct of cells)
load(stats_path, 'outputStats'); % per-subject MSClass etc.

OccurrenceNames = {'Occurrence_A', 'Occurrence_B', 'Occurrence_C', 'Occurrence_D', 'Occurrence_E', 'Occurrence_F', 'Occurrence_G'};

plotting = 0;

% Collect age vectors aligned to outputStats ordering

AgeMap = struct();
tmp = cellfun(@(s) s.age, France.BP_I.Depressed);
AgeMap.BP_I_Depressed = tmp(:);
tmp = cellfun(@(s) s.age, France.BP_I.Euthymic);
AgeMap.BP_I_Euthymic = tmp(:);
tmp = cellfun(@(s) s.age, France.BP_II.Depressed);
AgeMap.BP_II_Depressed = tmp(:);
tmp = cellfun(@(s) s.age, France.BP_II.Euthymic);
AgeMap.BP_II_Euthymic = tmp(:);
tmp = cellfun(@(s) s.age, France.HC);
AgeMap.HC = tmp(:);
tmp = cellfun(@(s) s.age, France.Siblings);
AgeMap.Siblings = tmp(:);

%% ----------------- Subject-level segments & transitions -------------- %%
% Out.<ERP>.<cond>.<group>{subj} fields:
% trans_counts (7x7; NaN diag), M_segments, N_transitions, C_segments (1x7), Occurrence_*.

Out = struct();
MicrostateLabels = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};

for e = 1:numel(ERPs)
    erp = ERPs{e};

    for c = 1:numel(Conds)
        condField = lower(Conds{c});

        for g = 1:numel(Groups)
            grp = Groups{g};
            subj_list = outputStats.(grp).(condField);
            nSubj = numel(subj_list);

            if isfield(AgeMap, grp) && numel(AgeMap.(grp)) ~= nSubj
                warning('%s: Age vector length (%d) != #subjects (%d).', grp, numel(AgeMap.(grp)), nSubj);
            end

            for sj = 1: nSubj
                S = subj_list(sj);

                if ~isfield(S,'MSClass') || isempty(S.MSClass)
                    continue;
                end

                segments = cell(1, size(S.MSClass, 2));
                num_transitions = zeros(1, size(S.MSClass, 2));
                Transition_matrix = zeros(7, 7);
                Transition_matrix(1: 8: 49) = NaN;

                % NOTE: assumes 'ERP.window{e}' exists in workspace

                for tri = 1:size(S.MSClass, 2)
                    aa = [S.MSClass(ERP.window{e}, tri); -999]; % sentinel
                    nonrepeats = find(diff(aa) ~= 0);
                    segments{tri} = aa(nonrepeats).';
                    segments{tri} = segments{tri}(segments{tri} ~= 0);
                    num_transitions(tri) = max(0, numel(segments{tri}) - 1);

                    for p = 1:(numel(segments{tri}) - 1)
                        i = segments{tri}(p);
                        j = segments{tri}(p + 1);

                        if i ~= j
                            Transition_matrix(i, j) = Transition_matrix(i, j) + 1;
                        end
                    end
                end

                M_segments = sum(cellfun(@numel, segments));
                pool = [segments{:}];
                C_segments = zeros(1, 7);

                for s = 1: 7
                    C_segments(s) = nnz(pool == s);
                end

                N_transitions = sum(num_transitions);

                Out.(erp).(condField).(grp){sj}.trans_counts = Transition_matrix;
                Out.(erp).(condField).(grp){sj}.M_segments = M_segments;
                Out.(erp).(condField).(grp){sj}.N_transitions = N_transitions;
                Out.(erp).(condField).(grp){sj}.C_segments = C_segments;
                for s = 1:7
                    Out.(erp).(condField).(grp){sj}.(OccurrenceNames{s}) = C_segments(s) / max(1, M_segments);
                end
                Out.(erp).(condField).(grp){sj}.transition_count_matrix = Transition_matrix / max(1, N_transitions);
            end
        end
    end
end

[srcIdx, tgtIdx] = find(~eye(7)); % 42 directed pairs
GLM_OPTS = statset('MaxIter', 5000, 'TolX', 1e-10, 'TolFun', 1e-10);

% Numerical guards (used below)

MIN_OBS = 3; % min pooled observed counts to fit a pair
MIN_EXP = 3; % min pooled expected (or anchored mean) to fit
LOGE_CLAMP = 12;  % clamp logE to [-12, +12] to stabilize IRLS

%% ============================================================= %%
%% Q1. Deviations from independence within each (ERP×Group×Condition)
%% ============================================================= %%

FDR_Q1_SCOPE = 'per_cell_42'; % options: 'per_cell_42' | 'pool_within_erp' | 'global_all'

TransStats_cellFDR = struct();

% For pooled Q1 options, collect p-values

poolQ1 = struct('p',[], 'map',[]);

for e = 1:numel(ERPs)
    erp = ERPs{e};

    for g = 1: numel(Groups)
        grp = Groups{g};

        % z-score age within group (for interpretability of intercept)

        Ages_raw = AgeMap.(grp)(:);
        muA = mean(Ages_raw, 'omitnan');
        sdA = std(Ages_raw, 'omitnan');

        if sdA <= eps
            AgesZ = zeros(size(Ages_raw));
        else
            AgesZ = (Ages_raw - muA)./sdA;
        end

        for c = 1: numel(Conds)
            condName = Conds{c};
            condField = lower(condName);
            cellOut = Out.(erp).(condField).(grp);

            IRR = NaN(7);
            beta0 = NaN(7);
            se = NaN(7);
            pval = NaN(7);
            ciLo = NaN(7);
            ciHi = NaN(7);
            nRows = zeros(7);
            ObsOverExp = NaN(7);

            % Long table for this cell

            RowsAll = table(); % y, logE, AgeZ, src, tgt

            if ~isempty(cellOut)
                for sj = 1:numel(cellOut)

                    if sj > numel(AgesZ)
                        continue;
                    end

                    S = cellOut{sj};

                    if isempty(S) || ~isfield(S,'trans_counts') || isempty(S.trans_counts)
                        continue;
                    end

                    if S.M_segments <= 0 || S.N_transitions <= 0
                        continue;
                    end

                    P = (S.C_segments(:) / S.M_segments);
                    E = S.N_transitions*(P * P.');
                    yvec = S.trans_counts(sub2ind([7, 7], srcIdx, tgtIdx));
                    evec = E(sub2ind([7, 7], srcIdx, tgtIdx));
                    keep = isfinite(evec) & evec > 0 & ~isnan(yvec);

                    if ~any(keep)
                        continue;
                    end

                    logE = log(evec(keep));
                    logE = max(min(logE, LOGE_CLAMP), -LOGE_CLAMP); % clamp

                    R = table(yvec(keep), logE, repmat(AgesZ(sj), nnz(keep), 1), srcIdx(keep), tgtIdx(keep), 'VariableNames', {'y', 'logE', 'AgeZ', 'src', 'tgt'});
                    RowsAll = [RowsAll; R];
                end
            end

            for k = 1: numel(srcIdx)
                s = srcIdx(k);
                t = tgtIdx(k);
                R = RowsAll(RowsAll.src==s & RowsAll.tgt == t, :);

                if isempty(R)
                    continue;
                end

                % screens

                if sum(R.y) < MIN_OBS || sum(exp(R.logE)) < MIN_EXP
                    continue;
                end

                nRows(s,t) = height(R);
                ObsOverExp(s,t) = sum(R.y) / sum(exp(R.logE));

                try
                    mdl = fitglm(R, 'y ~ 1 + AgeZ', 'Distribution', 'poisson', 'Link', 'log', 'Offset', R.logE, 'Options', GLM_OPTS);

                    b = mdl.Coefficients.Estimate;
                    seB = mdl.Coefficients.SE;
                    pB = mdl.Coefficients.pValue;
                    CI = coefCI(mdl, 0.05);

                    beta0(s, t) = b(1);
                    se(s, t) = seB(1);
                    pval(s, t) = pB(1);
                    IRR(s, t) = exp(b(1));
                    ciLo(s, t) = exp(CI(1, 1));
                    ciHi(s, t) = exp(CI(1, 2));

                    if ~strcmp(FDR_Q1_SCOPE, 'per_cell_42')
                        poolQ1.p(end + 1, 1)   = pval(s, t);
                        poolQ1.map(end + 1, :) = [e g c s t];
                    end
                catch
                    % leave NaNs
                end
            end

            % FDR for this cell (default)

            if strcmp(FDR_Q1_SCOPE, 'per_cell_42')
                [sigMask, qMat] = bh_fdr_mask(pval, 0.05);
            else
                sigMask = false(7);
                qMat = NaN(7);
            end

            dirMatrix = zeros(7);
            dirMatrix(sigMask & IRR > 1) =  1;
            dirMatrix(sigMask & IRR < 1) = -1;

            % Edge list

            [sSig, tSig] = find(sigMask);

            if ~isempty(sSig)
                irr_vec = IRR(sub2ind([7, 7], sSig, tSig));
                lo_vec = ciLo(sub2ind([7, 7], sSig, tSig));
                hi_vec = ciHi(sub2ind([7, 7], sSig, tSig));
                q_vec = qMat(sub2ind([7, 7], sSig, tSig));
                o2e_vec = ObsOverExp(sub2ind([7, 7], sSig, tSig));
                edges = table(sSig, tSig, irr_vec, lo_vec, hi_vec, o2e_vec, q_vec, 'VariableNames',{'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'ObsOverExp', 'qFDR'});
                edges.srcLabel = reshape(string(MicrostateLabels(edges.src)), [], 1);
                edges.tgtLabel = reshape(string(MicrostateLabels(edges.tgt)), [], 1);
                dlab = strings(numel(sSig), 1);
                dlab(irr_vec > 1) = "Above expected";
                dlab(irr_vec < 1) = "Below expected";
                edges.Direction = dlab;
            else
                edges = table([], [], [], [], [], [], [], 'VariableNames',{'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'ObsOverExp', 'qFDR'});
                edges.srcLabel = strings(0, 1);
                edges.tgtLabel = strings(0, 1);
                edges.Direction = strings(0, 1);
            end

            TransStats_cellFDR.(erp).(grp).(condName).IRR = IRR;
            TransStats_cellFDR.(erp).(grp).(condName).beta0 = beta0;
            TransStats_cellFDR.(erp).(grp).(condName).se = se;
            TransStats_cellFDR.(erp).(grp).(condName).p = pval;
            TransStats_cellFDR.(erp).(grp).(condName).q = qMat;
            TransStats_cellFDR.(erp).(grp).(condName).sig = sigMask;
            TransStats_cellFDR.(erp).(grp).(condName).dir = dirMatrix;
            TransStats_cellFDR.(erp).(grp).(condName).nRows = nRows;
            TransStats_cellFDR.(erp).(grp).(condName).ObsOverExp = ObsOverExp;
            TransStats_cellFDR.(erp).(grp).(condName).FDR_SCOPE = FDR_Q1_SCOPE;
            TransStats_cellFDR.(erp).(grp).(condName).edges = edges;
        end
    end
end

% Optional pooled BH for Q1

if ~strcmp(FDR_Q1_SCOPE, 'per_cell_42') && ~isempty(poolQ1.p)
    q_all = bh_adjust(poolQ1.p);

    % inject back

    for row = 1: size(poolQ1.map, 1)
        e = poolQ1.map(row, 1);
        g = poolQ1.map(row, 2);
        c = poolQ1.map(row, 3);
        s = poolQ1.map(row, 4);
        t = poolQ1.map(row, 5);
        erp = ERPs{e};
        grp = Groups{g};
        condName = Conds{c};
        qMat = TransStats_cellFDR.(erp).(grp).(condName).q;

        if all(size(qMat) == [1 1])
            qMat = NaN(7);
        end

        qMat(s, t) = q_all(row);
        TransStats_cellFDR.(erp).(grp).(condName).q = qMat;
    end

    % rebuild sig/dir/edges

    for e = 1: numel(ERPs)
        erp = ERPs{e};

        for g = 1: numel(Groups)
            grp = Groups{g};

            for c = 1:numel(Conds)
                condName = Conds{c};

                IRR = TransStats_cellFDR.(erp).(grp).(condName).IRR;
                qMat= TransStats_cellFDR.(erp).(grp).(condName).q;
                sig = qMat < 0.05;
                dir = zeros(7);
                dir(sig & IRR > 1) = 1;
                dir(sig & IRR < 1) = -1;
                TransStats_cellFDR.(erp).(grp).(condName).sig = sig;
                TransStats_cellFDR.(erp).(grp).(condName).dir = dir;

                [sSig, tSig] = find(sig);

                if ~isempty(sSig)
                    lo_vec = NaN(size(sSig));
                    hi_vec = NaN(size(sSig)); % CI optional
                    irr_vec = IRR(sub2ind([7, 7], sSig, tSig));
                    q_vec = qMat(sub2ind([7, 7], sSig, tSig));
                    o2e = TransStats_cellFDR.(erp).(grp).(condName).ObsOverExp(sub2ind([7,7], sSig, tSig));
                    edges = table(sSig,tSig, irr_vec, lo_vec, hi_vec, o2e, q_vec, 'VariableNames',{'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'ObsOverExp', 'qFDR'});
                    edges.srcLabel = reshape(string(MicrostateLabels(edges.src)), [], 1);
                    edges.tgtLabel = reshape(string(MicrostateLabels(edges.tgt)), [], 1);
                    dl = strings(numel(sSig), 1);
                    dl(irr_vec > 1) = "Above expected";
                    dl(irr_vec < 1) = "Below expected";
                    edges.Direction = dl;
                else
                    edges = table([], [], [], [], [], [], [], 'VariableNames',{'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'ObsOverExp', 'qFDR'});
                    edges.srcLabel = strings(0, 1);
                    edges.tgtLabel = strings(0,1);
                    edges.Direction = strings(0, 1);
                end
                TransStats_cellFDR.(erp).(grp).(condName).edges = edges;
            end
        end
    end
end

%% ======================================================= %%
% Q2. Anchored contrasts (Target vs HC at the same age)
%        FDR: within each Target group × ERP × Condition across the 42 pairs
%% ======================================================= %%

REF_GROUP = 'HC'; % <<< choose your anchor/reference group here
TARGET_GROUPS = {}; % {} => all groups except REF_GROUP; or subset, e.g. {'Siblings','BP_I_Euthymic'}

CompStats = buildAnchoredContrasts(Out, AgeMap, ERPs, Conds, Groups, 'RefGroup', REF_GROUP, 'TargetGroups', TARGET_GROUPS, 'Alpha', 0.05, 'MaxIter', 5000, 'MinObs', MIN_OBS, 'MinExp', MIN_EXP, 'LogEClamp', LOGE_CLAMP);
write_results_xlsx(save_path, ERPs, Conds, Groups, TransStats_cellFDR, CompStats, REF_GROUP);

%% ---------------------- GLOBAL WIDTH MAPPING ---------------------- %%
% Compute a global cap for |log(IRR)| across ALL Q1 and Q2 results

WidthDomain = computeGlobalWidthDomain(TransStats_cellFDR, CompStats, ERPs, Conds, Groups);

MIN_W = 0.3;
MAX_W = 11.0;
CURVE = 0.4;
NodeRadius = 0.2;
ArrowOffset = 0;
MicrostateLabels = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};
GreyColor = [0.90 0.90 0.90];
GreyLineWidth = 1.0;
ShaftGap = 0.16;
AboveColor = [0.85 0.10 0.10];
BelowColor = [0.10 0.35 0.95];

commonNV = {'NodeRadius',NodeRadius, 'CurveMagnitude', CURVE, 'WidthDomain', WidthDomain, 'MinWidth', MIN_W, 'MaxWidth', MAX_W, 'ArrowOffset', ArrowOffset, 'MicrostateLabels', MicrostateLabels, 'GreyColor', GreyColor, 'GreyLineWidth', GreyLineWidth, 'ShaftGap', ShaftGap, 'AboveColor', AboveColor, 'BelowColor', BelowColor};

%% ----------------------------- PLOTTING ------------------------------- %%

if plotting

    for e = 3: numel(ERPs)
        erp = ERPs{e};

        %%% Q1 (independence deviations) %%%%

        for g = 1: numel(Groups)
            grp = Groups{g};
            figName = sprintf('Q1 — Deviations from independence | %s — %s', erp, strrep(grp, '_', '  '));
            drawERPGridFigures(TransStats_cellFDR, erp, grp, Conds, MicrostateLabels, 'FigureName', figName, commonNV{:});
        end

        %%% Q2 (anchored contrasts) %%%%

        GroupContrastSlides(CompStats, erp, Conds, Groups, 'RefGroup', REF_GROUP, 'Alpha', 0.05, 'MicrostateLabels', MicrostateLabels, commonNV{:});
    end
end

function cap = sharedWidthCap_Q1Q2(TransStats, CompStats, erp, Conds, Groups, capDefault)

vals = [];

%%% Q1 IRRs %%%%

if isfield(TransStats, erp)
    for gi = 1: numel(Groups)
        grp = Groups{gi};

        if ~isfield(TransStats.(erp), grp)
            continue;
        end

        for ci = 1: numel(Conds)
            condName = Conds{ci};

            if ~isfield(TransStats.(erp).(grp), condName)
                continue;
            end

            node = TransStats.(erp).(grp).(condName);

            if isfield(node, 'IRR') && ~isempty(node.IRR)
                v = node.IRR(:);
                vals = [vals; v(isfinite(v) & v > 0)];
            end
        end
    end
end

%%% Q2 IRRs %%%%

if isfield(CompStats, erp)
    for ci = 1: numel(Conds)
        condName = Conds{ci};

        if ~isfield(CompStats.(erp), condName)
            continue;
        end

        S = CompStats.(erp).(condName);
        tg = fieldnames(S);

        for k = 1: numel(tg)
            node = S.(tg{k});

            if isfield(node, 'IRR') && ~isempty(node.IRR)
                v = node.IRR(:);
                vals = [vals; v(isfinite(v) & v > 0)];
            end
        end
    end
end

if isempty(vals)
    cap = capDefault;
else

    cap = min(capDefault, max(abs(log(vals))));

    if ~isfinite(cap) || cap <= 0
        cap = capDefault;
    end

end
end

function [sigMask, qMat] = bh_fdr_mask(pMat, alpha)

if nargin < 2
    alpha = 0.05;
end

v = pMat(:);
keep = ~isnan(v);
q = NaN(size(v));

if any(keep)
    q(keep) = bh_adjust(v(keep));
end

qMat = reshape(q, size(pMat));
sigMask = qMat < alpha;
end

function q = bh_adjust(p)

p = p(:);
[ps, ord] = sort(p, 'ascend');
m = numel(ps);
qtmp = ps .* m ./ (1: m)';
qtmp = flipud(cummin(flipud(qtmp)));
q = NaN(size(ps));
q(ord) = min(1, qtmp);
end

function CompStats = buildAnchoredContrasts(Out, AgeMap, ERPs, Conds, Groups, varargin)

% buildAnchoredContrasts (per-target 42-FDR version)
% Stage 1 (Ref-only, per pair): y ~ Poisson, log(mu) = log(E) + α0 + α_age*AgeS_ref
% Stage 2 (Ref+Target, per pair): y ~ Poisson, log(mu) = [log(E)+α0+α_age*AgeS_ref(row)] + δ*G (no intercept)
% FDR: BH within each ERP×Condition×TargetGroup across the 42 pairs.

p = inputParser;
addParameter(p, 'RefGroup', '', @(s)ischar(s)||isstring(s));
addParameter(p, 'TargetGroups', {}, @(x)iscell(x)||isstring(x));
addParameter(p, 'Alpha', 0.05);
addParameter(p, 'MaxIter', 2000);
addParameter(p, 'MinObs', 3);
addParameter(p, 'MinExp', 3);
addParameter(p, 'LogEClamp', 12);
parse(p, varargin{:});
opt = p.Results;

assert(~isempty(opt.RefGroup), 'Please provide ''RefGroup'', e.g., ''HC''.');
assert(any(strcmpi(Groups, opt.RefGroup)), 'RefGroup "%s" not found in Groups.', char(opt.RefGroup));

% Target list

if isempty(opt.TargetGroups)
    TargetList = Groups(~strcmpi(Groups, opt.RefGroup)); % all except reference
else
    tg = cellstr(opt.TargetGroups);
    TargetList = tg(~strcmpi(tg, opt.RefGroup) & ismember(tg, Groups));
end

assert(~isempty(TargetList), 'No valid TargetGroups to compare against RefGroup.');

[srcIdx, tgtIdx] = find(~eye(7)); % 42 directed pairs
GLM_OPTS = statset('MaxIter', opt.MaxIter, 'TolX', 1e-10, 'TolFun', 1e-10);
CompStats = struct();
MicrostateLabels = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};

for e = 1: numel(ERPs)
    erp = ERPs{e};

    for c = 1: numel(Conds)
        condName = Conds{c};
        condKey = lower(condName);

        RefCell = getfield_safe(Out, {erp, condKey, char(opt.RefGroup)});

        % Prebuild ref rows ONCE (shared across targets)

        Tref = table();
        if ~isempty(RefCell)
            agesRef = AgeMap.(char(opt.RefGroup));

            for sj = 1: numel(RefCell)
                S = RefCell{sj};

                if badSubject(S)
                    continue;
                end

                [yvec, evec] = obs_vs_exp(S, srcIdx, tgtIdx);
                keep = isfinite(evec) & evec > 0 & ~isnan(yvec);

                if ~any(keep)
                    continue;
                end

                logE = log(evec(keep));
                logE = max(min(logE, opt.LogEClamp), -opt.LogEClamp);

                R = table(yvec(keep), logE, repmat(agesRef(min(sj, end)), nnz(keep), 1), srcIdx(keep), tgtIdx(keep), 'VariableNames', {'y', 'logE', 'AgeRaw', 'src', 'tgt'});
                Tref = [Tref; R];
            end
        end

        for g = 1: numel(TargetList)
            tgtGrp = TargetList{g};
            Ttgt = table();
            TgtCell = getfield_safe(Out, {erp, condKey, tgtGrp});

            if ~isempty(TgtCell)
                agesTgt = AgeMap.(tgtGrp);

                for sj = 1: numel(TgtCell)
                    S = TgtCell{sj};

                    if badSubject(S)
                        continue;
                    end

                    [yvec, evec] = obs_vs_exp(S, srcIdx, tgtIdx);
                    keep = isfinite(evec) & evec > 0 & ~isnan(yvec);

                    if ~any(keep)
                        continue;
                    end

                    logE = log(evec(keep));
                    logE = max(min(logE, opt.LogEClamp), -opt.LogEClamp);

                    R = table(yvec(keep), logE, repmat(agesTgt(min(sj, end)), nnz(keep), 1), srcIdx(keep), tgtIdx(keep), 'VariableNames', {'y', 'logE', 'AgeRaw', 'src', 'tgt'});
                    Ttgt = [Ttgt; R];
                end
            end

            % Init per-target containers

            IRR = NaN(7);
            CIlo = NaN(7);
            CIhi = NaN(7);
            pMat = NaN(7);
            qMat = NaN(7);
            sig = false(7);

            if isempty(Tref) || isempty(Ttgt)
                store(erp, condName, tgtGrp, opt.RefGroup, IRR, CIlo, CIhi, pMat, qMat, sig, emptyEdges());
                continue;
            end

            % --- Fit per pair

            for k = 1: numel(srcIdx)
                s = srcIdx(k);
                t = tgtIdx(k);

                Rref = Tref(Tref.src == s & Tref.tgt == t, :);
                Rtgt = Ttgt(Ttgt.src == s & Ttgt.tgt == t, :);

                if isempty(Rref) || isempty(Rtgt)
                    continue;
                end

                % Stage 1: reference-only fit (age z-scored within ref)

                muR = mean(Rref.AgeRaw, 'omitnan');
                sdR = std(Rref.AgeRaw, 'omitnan');
                Rref_fit = Rref;

                if sdR <= eps
                    Rref_fit.AgeS = zeros(height(Rref_fit), 1);
                else
                    Rref_fit.AgeS = (Rref_fit.AgeRaw - muR)./sdR;
                end

                % screens for ref fit

                if sum(Rref_fit.y) < opt.MinObs || sum(exp(Rref_fit.logE)) < opt.MinExp
                    continue;
                end

                try
                    mdlRef = fitglm(Rref_fit, 'y ~ 1 + AgeS', 'Distribution', 'poisson', 'Link', 'log', 'Offset', Rref_fit.logE, 'Options', GLM_OPTS);

                    alpha0 = mdlRef.Coefficients.Estimate(1);
                    alphaAge = 0;

                    if any(strcmp(mdlRef.CoefficientNames, 'AgeS'))
                        alphaAge = mdlRef.Coefficients.Estimate(strcmp(mdlRef.CoefficientNames, 'AgeS'));
                    end
                catch
                    continue;
                end

                % Stage 2: combine ref+target & anchor offset at each row's age

                Rref2 = dropIfExists(Rref, {'AgeS', 'AgeS_ref', 'offset_ref', 'Gnum'});
                Rtgt2 = dropIfExists(Rtgt, {'AgeS', 'AgeS_ref', 'offset_ref', 'Gnum'});
                R2 = [Rref2; Rtgt2];

                if sdR <= eps
                    R2.AgeS_ref = zeros(height(R2), 1);
                else
                    R2.AgeS_ref = (R2.AgeRaw - muR)./sdR;
                end

                R2.offset_ref = R2.logE + alpha0 + alphaAge * R2.AgeS_ref;
                R2.Gnum = [zeros(height(Rref2), 1); ones(height(Rtgt2), 1)];

                % screens for combined fit

                if sum(R2.y) < opt.MinObs || sum(exp(R2.offset_ref)) < opt.MinExp
                    continue;
                end

                % quasi-separation guard

                y_t = sum(R2.y(R2.Gnum == 1));
                y_r = sum(R2.y(R2.Gnum == 0));

                if (y_t == 0 || y_r == 0)
                    continue;
                end

                try
                    mdl2 = fitglm(R2, 'y ~ -1 + Gnum', 'Distribution', 'poisson', 'Link', 'log', 'Offset', R2.offset_ref, 'Options', GLM_OPTS);

                    b = mdl2.Coefficients.Estimate(1);
                    p = mdl2.Coefficients.pValue(1);
                    CI = coefCI(mdl2, 0.05);

                    IRR(s, t)  = exp(b);
                    CIlo(s, t) = exp(CI(1, 1));
                    CIhi(s, t) = exp(CI(1, 2));
                    pMat(s, t) = p;

                catch
                    % leave NaNs
                end
            end

            % --- BH–FDR within this Target × ERP × Condition across 42 pairs

            offdiag = ~eye(7);
            p_lin = pMat(offdiag);
            keep = ~isnan(p_lin);

            if any(keep)
                q_all = bh_adjust(p_lin(keep));
                q_only = NaN(size(p_lin));
                q_only(keep) = q_all;

                qMat(offdiag) = q_only;
            end

            sig = qMat < opt.Alpha & ~isnan(qMat);

            % --- Edge list

            [sSig, tSig] = find(sig);

            if ~isempty(sSig)
                rr_vec = IRR(sub2ind([7, 7], sSig, tSig));
                lo_vec = CIlo(sub2ind([7, 7], sSig, tSig));
                hi_vec = CIhi(sub2ind([7, 7], sSig, tSig));
                q_vec  = qMat(sub2ind([7, 7], sSig, tSig));

                edges = table(sSig, tSig, rr_vec, lo_vec, hi_vec, q_vec, 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'qFDR'});
                edges.srcLabel = reshape(string(MicrostateLabels(edges.src)), [], 1);
                edges.tgtLabel = reshape(string(MicrostateLabels(edges.tgt)), [], 1);
                dlab = strings(numel(sSig), 1);
                dlab(rr_vec > 1) = "Higher vs " + string(opt.RefGroup);
                dlab(rr_vec < 1) = "Lower vs " + string(opt.RefGroup);
                edges.Direction = dlab;
            else
                edges = emptyEdges();
            end

            % --- Store

            store(erp, condName, tgtGrp, opt.RefGroup, IRR, CIlo, CIhi, pMat, qMat, sig, edges);
        end

        CompStats.(erp).(condName).FDR_SCOPE = 'per_group_condition_42pairs';
    end
end

    function tf = badSubject(S)
        tf = isempty(S) || ~isfield(S, 'trans_counts') || isempty(S.trans_counts) || S.M_segments <= 0 || S.N_transitions <= 0;
    end

    function [yvec, evec] = obs_vs_exp(S, srcIdx, tgtIdx)
        P = (S.C_segments(:) / S.M_segments); % 7×1
        E = S.N_transitions * (P * P.'); % 7×7 under independence
        yvec = S.trans_counts(sub2ind([7, 7], srcIdx, tgtIdx));
        evec = E(sub2ind([7, 7], srcIdx, tgtIdx));
    end

    function T = dropIfExists(T, names)
        keep = ismember(names, T.Properties.VariableNames);

        if any(keep)
            T = removevars(T, names(keep));
        end
    end

    function edges = emptyEdges()
        edges = table([], [], [], [], [], [], 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'qFDR'});
        edges.srcLabel = strings(0, 1);
        edges.tgtLabel = strings(0, 1);
        edges.Direction = strings(0, 1);
    end

    function S = getfield_safe(S, keys)
        for kk = 1: numel(keys)
            k = keys{kk};

            if ~isfield(S, k)
                S = [];
                return;
            end

            S = S.(k);
            if isempty(S)
                return;
            end

        end
    end

    function store(erp, condName, tgtGrp, refGrp, IRR, CIlo, CIhi, pMat, qMat, sig, edges)

        persistent CompStats_local
        if isempty(CompStats_local)
            CompStats_local = struct();
        end

        CompStats_local.(erp).(condName).(tgtGrp).IRR = IRR;
        CompStats_local.(erp).(condName).(tgtGrp).IRR_CIlo = CIlo;
        CompStats_local.(erp).(condName).(tgtGrp).IRR_CIhi = CIhi;
        CompStats_local.(erp).(condName).(tgtGrp).p = pMat;
        CompStats_local.(erp).(condName).(tgtGrp).q = qMat;
        CompStats_local.(erp).(condName).(tgtGrp).sig = sig;
        CompStats_local.(erp).(condName).(tgtGrp).edges = edges;
        CompStats_local.(erp).(condName).(tgtGrp).RefGroup = char(refGrp);
        assignin('caller', 'CompStats', CompStats_local);
    end

    function q = bh_adjust(p)
        p = p(:);
        [ps, ord] = sort(p, 'ascend');
        m = numel(ps);
        ranks = (1: m)';
        qtmp = (m ./ ranks) .* ps;
        qtmp = flipud(cummin(flipud(qtmp)));
        q = NaN(size(ps));
        q(ord) = min(1, qtmp);
    end
end

function WidthDomain = computeGlobalWidthDomain(TransStats_cellFDR, CompStats, ERPs, Conds, Groups)

% Returns [0, cap], where cap = max |log(IRR)| across:
% Q1: TransStats_cellFDR.(ERP).(Group).(Cond).IRR
% Q2: CompStats.(ERP).(Cond).(Target).IRR
% Off-diagonal pairs only; ignores NaNs/infs.

vals = [];

%%% Q1 %%%%

for e = 1: numel(ERPs)
    erp = ERPs{e};

    if ~isfield(TransStats_cellFDR, erp)
        continue;
    end

    for g = 1: numel(Groups)
        grp = Groups{g};

        if ~isfield(TransStats_cellFDR.(erp), grp)
            continue;
        end

        for c = 1: numel(Conds)
            cond = Conds{c};

            if ~isfield(TransStats_cellFDR.(erp).(grp), cond)
                continue;
            end

            node = TransStats_cellFDR.(erp).(grp).(cond);

            if ~isfield(node, 'IRR') || isempty(node.IRR)
                continue;
            end

            IRR = node.IRR;
            mask = ~eye(7) & isfinite(IRR);
            v = abs(log(IRR(mask)));
            vals = [vals; v(:)];
        end
    end
end

%%% Q2 %%%%

for e = 1: numel(ERPs)
    erp = ERPs{e};

    if ~isfield(CompStats, erp)
        continue;
    end

    for c = 1: numel(Conds)
        cond = Conds{c};

        if ~isfield(CompStats.(erp), cond)
            continue;
        end

        tgNames = fieldnames(CompStats.(erp).(cond));

        for k = 1: numel(tgNames)
            tgt = tgNames{k};
            node = CompStats.(erp).(cond).(tgt);

            if ~isfield(node, 'IRR') || isempty(node.IRR) || ~ismember(tgt, Groups)
                continue;
            end

            IRR = node.IRR;
            mask = ~eye(7) & isfinite(IRR);
            v = abs(log(IRR(mask)));
            vals = [vals; v(:)];
        end
    end
end

vals = vals(isfinite(vals) & vals >= 0);

if isempty(vals)
    cap = 1;
else
    cap = max(vals);
    if cap <= 0
        cap = 1;
    end
end

WidthDomain = [0 cap];
end
