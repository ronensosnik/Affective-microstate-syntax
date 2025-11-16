clear all
close all hidden
clc

% =======================================================================================================================================================================================
% A variables that will be used in all analyses: France
%
% The variable is a 1*1 struct with four fields, corresponding to group names: BP_I, BP_II, HC, Siblings. BP_I and BP_II, each, have two fields: Depressed and Euthymic (so, six groups)
% Each group field (e.g., HC, BP_I.Depressed, BP_II_Euthymic) is 1 * n cell, wherein n is the number of subjects in that group
%
% Each cell is 1*1 struct with 10 fields: number (double), age (double), gender (char), GAF (double), MADRS (double), YMRS (double), med (1*1 struct), Neutral (1*1 struct), Negative (1*1 struct), Positive (1*1 struct):
% med field is 1*1 struct with 9 fields (all double): AD, AP, MS, ANX, OTHER, poly_count, poly_ge2, drug_count, meds_missing)
% Neutral field is 1*1 struct with 6 fields: TooEarly (1 * 40 double), Correct_resp_RT (1 * 40 double), Incorrect_resp_RT (1 * 40 double), num_Correct_Incorrect_resp (1 * 40 double), num_no_resp (1 * 40 double), EGG_Data (1*1 with 43 fields)
% Negative field is 1*1 struct with 6 fields: TooEarly (1 * 40 double), Correct_resp_RT (1 * 40 double), Incorrect_resp_RT (1 * 40 double), num_Correct_Incorrect_resp (1 * 40 double), num_no_resp (1 * 40 double), EGG_Data (1*1 with 43 fields)
% Positive field is 1*1 struct with 6 fields: TooEarly (1 * 40 double), Correct_resp_RT (1 * 40 double), Incorrect_resp_RT (1 * 40 double), num_Correct_Incorrect_resp (1 * 40 double), num_no_resp (1 * 40 double), EGG_Data (1*1 with 43 fields)
% =======================================================================================================================================================================================

%% Microstate metrics

Pipeline_Sensitivity_tests = 0;
Meds_Sensitivity_tests = 0;

Find_microstates = 0;
Compute_and_save_metrics = 0;
Analyze_metrics_only_age_as_covariate = 0;
Analyze_transitions_only_age_as_covariate = 1;
Link_trial_RT = 0;

inputDir = 'D:\Papers\2025\Submitted\JAD (task microstates)\Data\';

if ~Pipeline_Sensitivity_tests
    saveDir = 'D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\';
    ClustPar_GFPPeaks = false; FitPar_PeakFit = 0;
else
    saveDir = 'D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Alternative pipelines\GFP_Peaks only\';
    ClustPar_GFPPeaks = true; FitPar_PeakFit = 1;
end

ERPs           = {'N200', 'P300', 'LPP'};
Microstates  = {'Microstate_A', 'Microstate_B', 'Microstate_C', 'Microstate_D', 'Microstate_E', 'Microstate_F', 'Microstate_G'};
Groups         = {'BP_I_Depressed', 'BP_I_Euthymic', 'BP_II_Depressed', 'BP_II_Euthymic', 'HC', 'Siblings'};
Conds          = {'Negative', 'Neutral', 'Positive'};
Occurrence  = {'Occurrence_A', 'Occurrence_B', 'Occurrence_C' ,'Occurrence_D' ,'Occurrence_E' ,'Occurrence_F', 'Occurrence_G'};

N200_electrodes_names = {'Fz', 'Cz', 'C1', 'C2', 'C3', 'C4', 'FC1', 'FC2', 'FC3', 'FC4'};
P300_electrodes_names = {'AFz', 'AF3', 'AF4', 'F3', 'Fz', 'F4'};
LPP_electrodes_names = {'CP1', 'CP2', 'CP3', 'CP4', 'P1', 'P2', 'P3', 'P4', 'PO3', 'PO4'};

ERP.window{1} = [181: 300];
ERP.window{2} = [301: 500];
ERP.window{3} = [501: 1000];

if Find_microstates

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
    Compare_grand_mean_Custo = 1;
    Backfit_and_quantify_temp_dynam = 1;

    %% Setting all parameters

    groupFolders = dir(inputDir);
    groupFolders = groupFolders([groupFolders.isdir]);
    groupFolders = groupFolders(~matches({groupFolders.name}, {'.', '..'}));
    nGroups = length(groupFolders);
    groupNames = cell(1, nGroups);
    dataDirs = cell(1, nGroups);

    % Set clustering parameters

    ClustPar.UseAAHC = false; % true = AAHC, false = kmeans
    ClustPar.MinClasses = 7; % minimum number of clusters to identify
    ClustPar.MaxClasses = 7; % maximum number of clusters to identify
    ClustPar.MaxMaps = inf; % maximum number of data samples to use to identify clusters
    ClustPar.GFPPeaks = ClustPar_GFPPeaks; % whether clustering should be limited to global field power peaks
    ClustPar.IgnorePolarity = true; % whether maps of inverted polarities should be considered part of the same cluster
    ClustPar.Normalize = true; % Set to false if using AAHC
    ClustPar.Allow_early_stop = true;

    ClustPar.Restarts = 100;

    % Set backfitting parameters

    FitPar.Classes = 7; % cluster solutions to use for backfitting
    FitPar.PeakFit = FitPar_PeakFit; % whether to backfit only on global field power peaks
    FitPar.lambda = 0.3; % smoothness penalty - ignored if FitPar.PeakFit = 1
    FitPar.b = 20; % smoothing window (ms) - ignored if FitPar.PeakFit = 1

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

            for j = 2: numel(dataDirs{i})

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

        if 0

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
        heatmap(sharedVarTable_plot, 'CellLabelColor', 'none', 'MissingDataColor', [1 1 1], 'GridVisible', 'off', 'ColorBarVisible' ,'off', 'Title', ['Comparison. Grand mean vs. Custo template. ' ERP_to_use_for_clustering.name]);
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
end

if Compute_and_save_metrics

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Within each window, labels were clipped to the window bounds.
    % For each microstate class we computed:
    % (i) coverage, the proportion of samples in the window carrying that label (segments spanning a boundary contributed by their overlap only);
    % (ii) occurrence, the number of onsets whose first sample fell within the window, divided by the window length to yield a rate (Hz);
    % (iii) mean duration, the average duration of segments whose onsets occurred within the window (segments truncated by the window boundary were excluded from the duration calculation).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    load('D:\Papers\2025\Submitted\JAD (task microstates)\Data\France_data.mat');
    load(fullfile(saveDir, ['outputStats.mat']), 'outputStats');

    % leaf template (the 5 fields at the bottom level)

    leaf = struct('num_occurrences', [], 'duration', [], 'coverage', [], 'number', [], 'age', [], 'meds', []);

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
                    Data.(e).(m).(g).(c) = repmat({leaf}, 1, 1);  % create the 6 fields here
                    Mean.(e).(m).(g).(c) = repmat({leaf}, 1, 1); % create the 6 fields here
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
                                y_ee(fff) = [];
                                run_starts_ee(fff) = [];
                                run_ends_ee(fff) = [];
                                z_ee(fff) = [];
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

                            % subject number, age and medications

                            switch Groups{gro}
                                case 'BP_I_Depressed'
                                    try
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_I.Depressed{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_I.Depressed{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds France.BP_I.Depressed{subj_num}.med];
                                    catch
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_I.Depressed{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_I.Depressed{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds France.BP_I.Depressed{subj_num}.med];
                                    end

                                case 'BP_II_Depressed'
                                    try
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_II.Depressed{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_II.Depressed{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_II.Depressed{subj_num}.med];
                                    catch
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_II.Depressed{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_II.Depressed{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds France.BP_II.Depressed{subj_num}.med];
                                    end

                                case 'BP_I_Euthymic'
                                    try
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_I.Euthymic{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_I.Euthymic{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds France.BP_I.Euthymic{subj_num}.med];
                                    catch
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_I.Euthymic{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_I.Euthymic{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds France.BP_I.Euthymic{subj_num}.med];
                                    end

                                case 'BP_II_Euthymic'
                                    try
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_II.Euthymic{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_II.Euthymic{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_II.Euthymic{subj_num}.med];
                                    catch
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.BP_II.Euthymic{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.BP_II.Euthymic{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds France.BP_II.Euthymic{subj_num}.med];
                                    end

                                case 'HC'
                                    try
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.HC{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.HC{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.HC{subj_num}.med];
                                    catch
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.HC{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.HC{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds France.HC{subj_num}.med];
                                    end

                                case 'Siblings'
                                    try
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.Siblings{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.Siblings{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.Siblings{subj_num}.med];
                                    catch
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age France.Siblings{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number France.Siblings{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds France.Siblings{subj_num}.med];
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
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = France.BP_I.Depressed{subj_num}.med;
                            case 'BP_II_Depressed'
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = France.BP_II.Depressed{subj_num}.number;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = France.BP_II.Depressed{subj_num}.age;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = France.BP_II.Depressed{subj_num}.med;
                            case 'BP_I_Euthymic'
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = France.BP_I.Euthymic{subj_num}.number;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = France.BP_I.Euthymic{subj_num}.age;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = France.BP_I.Euthymic{subj_num}.med;
                            case 'BP_II_Euthymic'
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = France.BP_II.Euthymic{subj_num}.number;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = France.BP_II.Euthymic{subj_num}.age;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = France.BP_II.Euthymic{subj_num}.med;
                            case 'HC'
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = France.HC{subj_num}.number;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = France.HC{subj_num}.age;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = France.HC{subj_num}.med;
                            case 'Siblings'
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.number = France.Siblings{subj_num}.number;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.age = France.Siblings{subj_num}.age;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conds{condd}){subj_num}.meds = France.Siblings{subj_num}.med;
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
            ALL_Mean.(eName).(mName).Meds                    = [];
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
                    meds = cellfun(@(s) s.meds, S, 'UniformOutput', false);
                    nums = cellfun(@(s) s.number, S);
                    covs = cellfun(@(s) s.coverage, S);
                    durs = cellfun(@(s) s.duration, S);
                    occs = cellfun(@(s) s.num_occurrences, S);
                    grps = cellfun(@(s) s.group, S, 'UniformOutput', false);
                    conds = cellfun(@(s) s.condition, S, 'UniformOutput', false);

                    % append (concatenate) to accumulators

                    ALL_Mean.(eName).(mName).Age                     = [ALL_Mean.(eName).(mName).Age,            ages];
                    ALL_Mean.(eName).(mName).Meds                   = [ALL_Mean.(eName).(mName).Meds,         meds];
                    ALL_Mean.(eName).(mName).Number               = [ALL_Mean.(eName).(mName).Number,      nums];
                    ALL_Mean.(eName).(mName).Coverage            = [ALL_Mean.(eName).(mName).Coverage,   covs];
                    ALL_Mean.(eName).(mName).Duration              = [ALL_Mean.(eName).(mName).Duration,     durs];
                    ALL_Mean.(eName).(mName).Num_occurrence = [ALL_Mean.(eName).(mName).Num_occurrence, occs];
                    ALL_Mean.(eName).(mName).Group                  = [ALL_Mean.(eName).(mName).Group,         grps];
                    ALL_Mean.(eName).(mName).Condition             = [ALL_Mean.(eName).(mName).Condition,    conds];
                end
            end
        end
    end

    if ~Pipeline_Sensitivity_tests
        save('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Metrics_ALL.mat', 'ALL_Mean');
    else
        save('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Alternative pipelines\GFP_Peaks only\Metrics_ALL', 'ALL_Mean')
    end
end

if Analyze_metrics_only_age_as_covariate

    if ~Pipeline_Sensitivity_tests
        load('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Metrics_ALL.mat');
    else
        load('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Alternative pipelines\GFP_Peaks only\Metrics_ALL.mat');
    end

    % For each of the three measures, we run 21 sets of models (7 microstates * 3 ERPs)

    if ~Pipeline_Sensitivity_tests
        Out_Duration_base = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Duration', 0.05);
        Out_Occurrence_base = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Occurrence', 0.05);
        Out_Coverage_base = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Coverage', 0.05);
    else
        Out_Duration_sens_peaks = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Duration', 0.05);
        Out_Occurrence_sens_peaks = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Occurrence', 0.05);
        Out_Coverage_sens_peaks = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Coverage', 0.05);
    end

    % Report

    if ~Pipeline_Sensitivity_tests
        print_summary_measure(Out_Duration_base, 0.05);
        print_summary_measure(Out_Occurrence_base, 0.05);
        print_summary_measure(Out_Coverage_base, 0.05);
    else
        print_summary_measure(Out_Duration_sens_peaks, 0.05);
        print_summary_measure(Out_Occurrence_sens_peaks, 0.05);
        print_summary_measure(Out_Coverage_sens_peaks, 0.05);
    end

    % Generate plots for any FDR-significant effects

    if ~Pipeline_Sensitivity_tests
        plot_emm_significant(Out_Duration_base,      'Duration',      0.05, ALL_Mean);
        plot_emm_significant(Out_Occurrence_base, 'Occurrence', 0.05, ALL_Mean);
        plot_emm_significant(Out_Coverage_base,    'Coverage',    0.05, ALL_Mean);
    else
        plot_emm_significant(Out_Duration_sens_peaks,      'Duration',      0.05, ALL_Mean);
        plot_emm_significant(Out_Occurrence_sens_peaks, 'Occurrence', 0.05, ALL_Mean);
        plot_emm_significant(Out_Coverage_sens_peaks,    'Coverage',    0.05, ALL_Mean);

    end

    if ~Pipeline_Sensitivity_tests
        save(fullfile(saveDir, ['Out_Duration_base.mat']), 'Out_Duration_base');
        save(fullfile(saveDir, ['Out_Occurrence_base.mat']), 'Out_Occurrence_base');
        save(fullfile(saveDir, ['Out_Coverage_base.mat']), 'Out_Coverage_base');

        export_pairwise_S5(Out_Duration_base, 'Duration', 'baseline', 'S5_pairs_baseline_Dur.csv');
        export_pairwise_S5(Out_Occurrence_base, 'Occurrence', 'baseline', 'S5_pairs_baseline_Occ.csv');
        export_pairwise_S5(Out_Coverage_base, 'Coverage', 'baseline', 'S5_pairs_baseline_Cov.csv');
    else
        save(fullfile(saveDir, ['Out_Duration_sens_peaks.mat']), 'Out_Duration_sens_peaks');
        save(fullfile(saveDir, ['Out_Occurrence_sens_peaks.mat']), 'Out_Occurrence_sens_peaks');
        save(fullfile(saveDir, ['Out_Coverage_sens_peaks.mat']), 'Out_Coverage_sens_peaks');

        export_pairwise_S5(Out_Duration_sens_peaks, 'Duration',  'peaks_only', 'S5_pairs_peaks_Dur.csv');
        export_pairwise_S5(Out_Occurrence_sens_peaks, 'Occurrence',  'peaks_only', 'S5_pairs_peaks_Occ.csv');
        export_pairwise_S5(Out_Coverage_sens_peaks, 'Coverage',  'peaks_only', 'S5_pairs_peaks_Cov.csv');
    end

end

if Analyze_transitions_only_age_as_covariate

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
    %%% ====================================================================== %%

    %% ------------------------------ I/O ---------------------------------- %%

    data_path = 'D:\Papers\2025\Submitted\JAD (task microstates)\Data\France_data.mat';
    stats_path = [saveDir '\outputStats.mat'];

    load(data_path, 'France'); % ages & metadata (struct of cells)
    load(stats_path, 'outputStats'); % per-subject MSClass etc.

    OccurrenceNames = {'Occurrence_A', 'Occurrence_B', 'Occurrence_C', 'Occurrence_D', 'Occurrence_E', 'Occurrence_F', 'Occurrence_G'};

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

    for e = 1: numel(ERPs)
        erp = ERPs{e};

        for c = 1: numel(Conds)
            condField = lower(Conds{c});

            for g = 1: numel(Groups)
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

                    for tri = 1: size(S.MSClass, 2)
                        aa = [S.MSClass(ERP.window{e}, tri); -999]; % sentinel
                        nonrepeats = find(diff(aa) ~= 0);
                        segments{tri} = aa(nonrepeats).';
                        segments{tri} = segments{tri}(segments{tri} ~= 0); % number of labels in a trial
                        num_transitions(tri) = max(0, numel(segments{tri}) - 1);

                        for p = 1:(numel(segments{tri}) - 1)
                            i = segments{tri}(p);
                            j = segments{tri}(p + 1);

                            if i ~= j
                                Transition_matrix(i, j) = Transition_matrix(i, j) + 1;
                            end
                        end
                    end

                    M_segments = sum(cellfun(@numel, segments)); % total number of labels (in all trials)
                    pool = [segments{:}]; % a list of all labels
                    C_segments = zeros(1, 7);

                    for s = 1: 7
                        C_segments(s) = nnz(pool == s); % total number of each lf the 7 labels (in all trials)
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

    LOGE_CLAMP = 12; % clamp logE to [-12, +12] to stabilize IRLS
    MIN_OBS = 3; % min pooled observed counts to fit a pair
    MIN_EXP = 3; % min pooled expected (or anchored mean) to fit

    %% =========================================================== %%
    %% Q1. Deviations from independence within each (ERP×Group×Condition)
    %% =========================================================== %%

    FDR_Q1_SCOPE = 'per_cell_42'; % options: 'per_cell_42' | 'pool_within_erp' | 'global_all'

    TransStats_cellFDR = struct();

    % For pooled Q1 options, collect p-values

    poolQ1 = struct('p', [], 'map',[]);

    for e = 1: numel(ERPs)
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
                    for sj = 1: numel(cellOut)

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

                        %  The below code does not take in consideration and omits self transitions (e.g., A -> A) so the expected count is smaller and it can lead to more  transitions than deviate from independence

                        if 0
                            E = S.N_transitions*(P * P.');
                            yvec = S.trans_counts(sub2ind([7, 7], srcIdx, tgtIdx));
                            evec = E(sub2ind([7, 7], srcIdx, tgtIdx));
                            keep = isfinite(evec) & evec > 0 & ~isnan(yvec);
                        end

                        %  The below code DOES take in consideration and omits self transitions (e.g., A -> A)

                        if 1
                            Ssq = sum(P.^2);
                            den = 1 - Ssq;
                            kK = numel(S.C_segments);

                            gamma = 1 / den;

                            E = gamma * S.N_transitions * (P * P.');
                            E(1: kK + 1: kK  *kK) = NaN;

                            yvec = S.trans_counts(sub2ind([kK, kK], srcIdx, tgtIdx));
                            evec = E(sub2ind([kK, kK], srcIdx, tgtIdx));
                            keep = isfinite(evec) & evec > 0 & ~isnan(yvec);
                        end

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

                    nRows(s, t) = height(R);
                    ObsOverExp(s, t) = sum(R.y) / sum(exp(R.logE));

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
                            poolQ1.p(end + 1, 1) = pval(s, t);
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
                dirMatrix(sigMask & IRR > 1) = 1;
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

    %% ======================================================= %%
    % Q2. Anchored contrasts (Target vs HC at the same age)
    %        FDR: within each Target group × ERP × Condition across the 42 pairs
    %% ======================================================= %%

    REF_GROUP = 'HC'; % <<< choose your anchor/reference group here
    TARGET_GROUPS = {}; % {} => all groups except REF_GROUP; or subset, e.g. {'Siblings','BP_I_Euthymic'}

    CompStats = buildAnchoredContrasts(Out, AgeMap, ERPs, Conds, Groups, 'RefGroup', REF_GROUP, 'TargetGroups', TARGET_GROUPS, 'Alpha', 0.05, 'MaxIter', 5000, 'MinObs', MIN_OBS, 'MinExp', MIN_EXP, 'LogEClamp', LOGE_CLAMP);
    write_results_xlsx(saveDir, ERPs, Conds, Groups, TransStats_cellFDR, CompStats, REF_GROUP);

    %% ----------------------------- PLOTTING ------------------------------- %%

    plotting = 1;

    if plotting

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

        for e = 1: numel(ERPs)

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
end

if Meds_Sensitivity_tests

    %% ====================================== %
    %  Static measures (Duration / Occurrence / Coverage)
    %    =======================================%

    data_path = 'D:\Papers\2025\Submitted\JAD (task microstates)\Data\France_data.mat';
    stats_path = 'D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\outputStats.mat';
    save_path = 'D:\Papers\2025\Submitted\JAD (task microstates)\Analysis';

    load('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Metrics_ALL.mat');

    % For each of the three measures, we run 21 sets of models (7 microstates * 3 ERPs)

    Out_Duration_Base = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Duration', 0.05);
    Out_Occurrence_Base = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Occurrence', 0.05);
    Out_Coverage_Base = lme_posthoc_21FDR_transformsBT(ALL_Mean, 'Coverage', 0.05);

    Out_Duration_med = lme_posthoc_21FDR_transformsBT_sens(ALL_Mean, 'Duration', 0.05);
    Out_Occurrence_med = lme_posthoc_21FDR_transformsBT_sens(ALL_Mean, 'Occurrence', 0.05);
    Out_Coverage_med = lme_posthoc_21FDR_transformsBT_sens(ALL_Mean, 'Coverage', 0.05);

    Out_Duration_med_poly = lme_posthoc_21FDR_transformsBT_sens(ALL_Mean, 'Duration', 0.05, 'poly');
    Out_Occurrence_med_poly = lme_posthoc_21FDR_transformsBT_sens(ALL_Mean, 'Occurrence', 0.05, 'poly');
    Out_Coverage_med_poly = lme_posthoc_21FDR_transformsBT_sens(ALL_Mean, 'Coverage', 0.05, 'poly');

    % Report sensitivity results

    Duration_Baseline_vs_meds = summarize_metric_meds_robustness(Out_Duration_Base, Out_Duration_med, 'Duration', 0.05)
    Occurrence_Baseline_vs_meds = summarize_metric_meds_robustness(Out_Occurrence_Base, Out_Occurrence_med, 'Occurrence', 0.05)
    Coverage_Baseline_vs_meds = summarize_metric_meds_robustness(Out_Coverage_Base, Out_Coverage_med, 'Coverage', 0.05)

    Duration_Baseline_vs_meds_poly =  summarize_metric_meds_robustness(Out_Duration_Base, Out_Duration_med_poly, 'Duration', 0.05)
    Occurrence_Baseline_vs_meds_poly = summarize_metric_meds_robustness(Out_Occurrence_Base, Out_Occurrence_med_poly, 'Occurrence', 0.05)
    Coverage_Baseline_vs_meds_poly = summarize_metric_meds_robustness(Out_Coverage_Base, Out_Coverage_med_poly, 'Coverage', 0.05)

    save('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Duration_Baseline_vs_meds.mat', 'Duration_Baseline_vs_meds');
    save('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Occurrence_Baseline_vs_meds.mat', 'Occurrence_Baseline_vs_meds');
    save('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Coverage_Baseline_vs_meds.mat', 'Coverage_Baseline_vs_meds');

    save('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Duration_Baseline_vs_meds_poly.mat', 'Duration_Baseline_vs_meds_poly');
    save('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Occurrence_Baseline_vs_meds_poly.mat', 'Occurrence_Baseline_vs_meds_poly');
    save('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Coverage_Baseline_vs_meds_poly.mat', 'Coverage_Baseline_vs_meds_poly');

    convert_meds_robustness_to_excel; % auto-discovers the 6 files and writes Metrics_Meds_Robustness.xlsx

    load(data_path, 'France'); % ages & metadata (struct of cells)
    load(stats_path, 'outputStats'); % per-subject MSClass etc.

    %% ----------------------- Fixed settings / labels --------------------- %%

    ERPs = {'N200', 'P300', 'LPP'};
    Groups = {'BP_I_Depressed', 'BP_I_Euthymic', 'BP_II_Depressed', 'BP_II_Euthymic', 'HC', 'Siblings'};
    Conds = {'Neutral', 'Negative', 'Positive'};
    MicrostateLabels = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};
    OccurrenceNames = {'Occurrence_A', 'Occurrence_B', 'Occurrence_C', 'Occurrence_D', 'Occurrence_E', 'Occurrence_F', 'Occurrence_G'};

    % FDR scope: BH across 42 directed edges per cell (7x7 off-diagonal)

    FDR_Q1_SCOPE = 'per_cell_42';

    % GLM & safety clamps

    GLM_OPTS = statset('MaxIter', 5000, 'TolX', 1e-10, 'TolFun', 1e-10);
    LOGE_CLAMP = 30; % clamp for log(E)
    MIN_OBS = 3; % min observed transition count to fit a model for an edge
    MIN_EXP = 3; % min expected transitions for an edge

    alphaFDR = 0.05;

    %% ----------------- Build AgeMap aligned to outputStats --------------- %%

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
    % Out.<ERP>.<cond>.<group>{subj}:
    %   trans_counts (7x7; NaN diag), M_segments, N_transitions,
    %   C_segments (1x7), Occurrence_*, transition_count_matrix (normalized)

    Out = struct();

    for e = 1: numel(ERPs)
        erp = ERPs{e};
        condFields = lower(Conds);

        for c = 1: numel(Conds)
            condField = condFields{c};

            for g = 1: numel(Groups)
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

                    winIdx = get_erp_window_idx(e, S); % uses ERP.window{e} if present, else all rows
                    segments = cell(1, size(S.MSClass, 2));
                    num_transitions = zeros(1, size(S.MSClass, 2));

                    Transition_matrix = zeros(7, 7);
                    Transition_matrix(1:8:49) = NaN;

                    for tri = 1: size(S.MSClass, 2)
                        aa = [S.MSClass(winIdx, tri); -999]; % sentinel
                        nonrepeats = find(diff(aa) ~= 0);
                        segments{tri} = aa(nonrepeats).';
                        segments{tri} = segments{tri}(segments{tri} ~= 0);
                        num_transitions(tri) = max(0, numel(segments{tri}) - 1);

                        for p = 1: (numel(segments{tri}) - 1)
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

                    for s = 1: 7
                        Out.(erp).(condField).(grp){sj}.(OccurrenceNames{s}) = C_segments(s) / max(1, M_segments);
                    end
                    Out.(erp).(condField).(grp){sj}.transition_count_matrix = Transition_matrix / max(1, N_transitions);
                end
            end
        end
    end

    [srcIdx, tgtIdx] = find(~eye(7)); % 42 directed pairs

    %% ===================================================
    %% Q1 BASELINE (Age-only), BH–FDR within each cell (42 edges)
    %% ===================================================
    TransStats_cellFDR = q1_core(Out, AgeMap, ERPs, Groups, Conds, FDR_Q1_SCOPE, GLM_OPTS, LOGE_CLAMP, MIN_OBS, MIN_EXP, MicrostateLabels, 'age_only', alphaFDR);

    save(fullfile(save_path,'TransStats_Q1_Baseline_FDR.mat'), 'TransStats_cellFDR');

    %% ===================================================
    %% Q1 with Medication Covariates (two sensitivity specifications)
    %% ===================================================

    % Build MedMap (aligned to group subject order)

    MedMap = build_medmap(France, Groups);

    % (A) class dummies (AD/AP/MS/ANX/OTHER)

    TransStats_cellFDR_classes = q1_core(Out, AgeMap, ERPs, Groups, Conds, FDR_Q1_SCOPE, GLM_OPTS, LOGE_CLAMP, MIN_OBS, MIN_EXP, MicrostateLabels, 'classes', alphaFDR, MedMap);
    save('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Including_meds_as_covariates\Transitions\Q1\TransStats_Q1_MedsClasses_FDR.mat', 'TransStats_cellFDR_classes');

    % (B) polypharmacy only (poly_count)

    TransStats_cellFDR_poly = q1_core(Out, AgeMap, ERPs, Groups, Conds, FDR_Q1_SCOPE, GLM_OPTS, LOGE_CLAMP, MIN_OBS, MIN_EXP, MicrostateLabels, 'poly', alphaFDR, MedMap);
    save('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Including_meds_as_covariates\Transitions\Q1\TransStats_Q1_MedsPoly_FDR.mat', 'TransStats_cellFDR_poly');

    [overall_A, byERP_A] = summarize_q1_robustness(TransStats_cellFDR, TransStats_cellFDR_classes, ERPs, Groups, Conds);
    [overall_B, byERP_B] = summarize_q1_robustness(TransStats_cellFDR, TransStats_cellFDR_poly, ERPs, Groups, Conds);

    fprintf('\nQ1 robustness — classes (AD/AP/MS/ANX/OTHER):\n'); disp(overall_A);
    fprintf('\nBy ERP:\n');
    disp(byERP_A);

    fprintf('\nQ1 robustness — poly_count only:\n');
    disp(overall_B);
    fprintf('\nBy ERP:\n');
    disp(byERP_B);

    xlsx_q1 = fullfile('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Including_meds_as_covariates\Transitions\Q1\', 'Q1_Meds_Robustness.xlsx');

    if exist(xlsx_q1, 'file') == 2
        delete(xlsx_q1);
    end

    writetable(struct2table(overall_A), xlsx_q1, 'Sheet','classes_overall');
    writetable(byERP_A, xlsx_q1, 'Sheet', 'classes_by_ERP');
    writetable(struct2table(overall_B), xlsx_q1, 'Sheet', 'poly_overall');
    writetable(byERP_B, xlsx_q1, 'Sheet', 'poly_by_ERP');

    fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', fullfile('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Including_meds_as_covariates\Transitions\Q1\TransStats_Q1_Baseline_FDR.mat'), fullfile(save_path,'TransStats_Q1_MedsClasses_FDR.mat'), fullfile('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Including_meds_as_covariates\Transitions\Q1\\ransStats_Q1_MedsPoly_FDR.mat'), xlsx_q1);

    %% ===================== %%
    %% Q2 BASELINE (Age-only)
    %% ===================== %%

    REF_GROUP = 'HC'; % anchor/reference
    TARGET_GROUPS = {}; % {} => all except REF; or subset, e.g., {'Siblings', 'BP_I_Euthymic'}
    alphaFDR = 0.05;
    MAXITER = 5000;

    % ---------- Baseline Q2 (age-only; uses your existing function) ----------

    CompStats_Q2 = buildAnchoredContrasts(Out, AgeMap, ERPs, Conds, Groups, 'RefGroup', REF_GROUP, 'TargetGroups', TARGET_GROUPS, 'Alpha', alphaFDR, 'MaxIter', MAXITER, 'MinObs', MIN_OBS, 'MinExp', MIN_EXP, 'LogEClamp', LOGE_CLAMP);

    save('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Including_meds_as_covariates\Transitions\Q2\CompStats_Q2_Baseline_FDR.mat', 'CompStats_Q2');

    % ---------- Build MedMap aligned to group subject order ----------

    MedMap = build_medmap_for_q2(France, Groups);

    % ---------- Q2 with meds (A) psychotropic classes ----------

    CompStats_Q2_classes = buildAnchoredContrasts_meds(Out, AgeMap, MedMap, ERPs, Conds, Groups, 'RefGroup', REF_GROUP, 'TargetGroups', TARGET_GROUPS, 'Spec', 'classes', 'Alpha', alphaFDR, 'MaxIter', MAXITER, 'MinObs', MIN_OBS, 'MinExp', MIN_EXP, 'LogEClamp', LOGE_CLAMP);

    save('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Including_meds_as_covariates\Transitions\Q2\CompStats_Q2_MedsClasses_FDR.mat', 'CompStats_Q2_classes');

    % ---------- Q2 with meds (B) polypharmacy only ----------

    CompStats_Q2_poly = buildAnchoredContrasts_meds(Out, AgeMap, MedMap, ERPs, Conds, Groups, 'RefGroup', REF_GROUP, 'TargetGroups', TARGET_GROUPS, 'Spec', 'poly', 'Alpha', alphaFDR, 'MaxIter', MAXITER, 'MinObs', MIN_OBS, 'MinExp', MIN_EXP, 'LogEClamp', LOGE_CLAMP);

    save('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Including_meds_as_covariates\Transitions\Q2\CompStats_Q2_MedsPoly_FDR.mat', 'CompStats_Q2_poly');

    % ---------- Summarize robustness vs baseline and write Excel ----------

    [overall_A, byERP_A] = summarize_q2_robustness(CompStats_Q2, CompStats_Q2_classes, ERPs, Conds, Groups);
    [overall_B, byERP_B] = summarize_q2_robustness(CompStats_Q2, CompStats_Q2_poly, ERPs, Conds, Groups);

    fprintf('\nQ2 robustness — classes (AD/AP/MS/ANX/OTHER):\n');
    disp(overall_A);

    fprintf('\nBy ERP:\n');
    disp(byERP_A);

    fprintf('\nQ2 robustness — poly_count only:\n');
    disp(overall_B);

    fprintf('\nBy ERP:\n');
    disp(byERP_B);

    xlsx_q2 = fullfile(save_path, 'Q2_Meds_Robustness.xlsx');

    if exist(xlsx_q2, 'file') == 2
        delete(xlsx_q2);
    end

    writetable(struct2table(overall_A), xlsx_q2, 'Sheet', 'classes_overall');
    writetable(byERP_A, xlsx_q2, 'Sheet', 'classes_by_ERP');
    writetable(struct2table(overall_B), xlsx_q2, 'Sheet','poly_overall');
    writetable(byERP_B, xlsx_q2, 'Sheet','poly_by_ERP');

    fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', fullfile('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Including_meds_as_covariates\Transitions\Q2\CompStats_Q2_Baseline_FDR.mat'), fullfile('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Including_meds_as_covariates\Transitions\Q2\CompStats_Q2_MedsClasses_FDR.mat'), fullfile('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Including_meds_as_covariates\Transitions\Q2\CompStats_Q2_MedsPoly_FDR.mat'), xlsx_q2);

    %% =========================================
    %  Q3 — Behavior–syntax link (with meds)
    %  Confirmatory family (Anchor-2 slopes)
    %  Saves: rt_trial_confirmatory2_meds_[classes|poly].csv
    %% =========================================

    % Reuse demographics and per-trial RT + syntax features

    load(data_path, 'France'); % ages & meds live here
    S = load(stats_path);
    outputStats = S.outputStats; % MSClass etc.

    % ERP windows (match your confirmatory family)

    windows.N200 = [180 300];
    windows.P300 = [300 500];
    windows.LPP = [500 1000];

    % HC-derived anchors and LPP backbone

    anchors = findAnchorsBackbone(outputStats, windows);

    % Map subject numbers to outputStats order

    subj_id_map = buildSubjectIdMap(France);

    % Trial-level RT table

    rt_trials = buildRTTrialsTable(France); % subject, condition, trial, rt_ms, log_rt, age

    % Trial-level syntax predictors (rates within ERP windows)

    trial_feats = computeTrialPredictors(outputStats, anchors, windows, subj_id_map);

    % Join and prepare

    trial_data = innerjoin(rt_trials, trial_feats, 'Keys', {'subject', 'condition', 'trial'});
    trial_data.condition = categorical(trial_data.condition, {'Neutral', 'Negative', 'Positive'});

    % Standardize predictors globally

    pred_cols = {'N200_a1_rate', 'N200_a2_rate', 'P300_a1_rate', 'P300_a2_rate', 'LPP_backbone_rate'};

    for k = 1: numel(pred_cols)
        zname = ['z_' pred_cols{k}];
        mu = mean(trial_data.(pred_cols{k}), 'omitnan');
        sd = std(trial_data.(pred_cols{k}), 'omitnan');

        if sd==0 || isnan(sd), trial_data.(zname) = zeros(height(trial_data), 1);
        else
            trial_data.(zname) = (trial_data.(pred_cols{k}) - mu) ./ sd;
        end
    end

    trial_data.trial_z = (trial_data.trial - mean(trial_data.trial)) ./ std(trial_data.trial);

    % Within subject×condition demeaning for the three predictors we test

    pred_z = {'z_N200_a2_rate', 'z_P300_a2_rate', 'z_LPP_backbone_rate'};
    trial_data = demeanWithinSubjCond(trial_data, pred_z);

    % Keep complete rows for confirmatory family

    use = isfinite(trial_data.log_rt) & isfinite(double(trial_data.age)) & isfinite(trial_data.trial_z);
    need = {'z_N200_a2_rate', 'z_P300_a2_rate', 'z_LPP_backbone_rate'};

    for k = 1: numel(need)
        use = use & isfinite(trial_data.(['w_' need{k}])) & isfinite(trial_data.(['m_' need{k}]));
    end

    trialW = trial_data(use, :);
    trialW.c_val = double(trialW.condition ~= 'Neutral'); % Neutral=0; Neg/Pos=1

    % ---------- Two medication specifications: (A) five classes; (B) poly only ----------

    for spec = ["classes", "poly"]

        % Build subject-level meds table and left-join (fill missing with 0)

        medsT = build_meds_subject_table(France, char(spec)); % columns: subject + meds
        trialWm = outerjoin(trialW, medsT, 'Keys', 'subject', 'MergeKeys', true, 'Type', 'left');

        varNames = trialWm.Properties.VariableNames; % cellstr or string array

        for i = 1: numel(varNames)
            vn = varNames{i}; % <- extract to char

            if isnumeric(trialWm.(vn)) || islogical(trialWm.(vn))
                trialWm.(vn) = double(trialWm.(vn));
            end
        end

        medVars = setdiff(medsT.Properties.VariableNames, {'subject'});

        % Fit the two confirmatory models (Anchor-2) + the LPP control, now controlling for meds

        Rn200 = fit_one_anchor2_meds(trialWm, 'z_N200_a2_rate', 'N200 Anchor-2 (valenced vs neutral)', medVars);
        Rp300 = fit_one_anchor2_meds(trialWm, 'z_P300_a2_rate', 'P300 Anchor-2 (valenced vs neutral)', medVars);
        Rctl = fit_one_anchor2_meds(trialWm, 'z_LPP_backbone_rate', 'LPP backbone (control; not counted)', medVars);

        confirm2 = [Rn200; Rp300];
        confirm2.p_Holm = holm_stepdown(confirm2.pValue); % two tests

        % Ensure Rctl has the same vars as confirm2 (esp. p_Holm)

        vars = confirm2.Properties.VariableNames;

        % Add any missing variables to Rctl with appropriate defaults

        for k = 1: numel(vars)
            vn = vars{k};

            if ~ismember(vn, Rctl.Properties.VariableNames)

                % create a default column that matches the class used in confirm2

                ex = confirm2.(vn);

                if isnumeric(ex)
                    Rctl.(vn) = NaN(height(Rctl), 1);
                elseif isstring(ex)
                    Rctl.(vn) = strings(height(Rctl), 1);
                    Rctl.(vn)(:) = missing;
                elseif islogical(ex)
                    Rctl.(vn) = false(height(Rctl), 1);   % logical has no NaN; false is fine
                elseif iscategorical(ex)
                    Rctl.(vn) = categorical( strings(height(Rctl), 1), categories(ex));
                else
                    % fallback: replicate first element's type/size
                    Rctl.(vn) = repmat(ex(1), height(Rctl), 1);
                end
            end
        end

        % If confirm2 is missing any vars that Rctl has (unlikely), drop them for alignment

        extraInRctl = setdiff(Rctl.Properties.VariableNames, vars, 'stable');
        if ~isempty(extraInRctl)
            Rctl = removevars(Rctl, extraInRctl);
        end

        % Reorder columns to match confirm2

        Rctl = Rctl(:, vars);

        outT = [confirm2; Rctl];
        out_csv = fullfile('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Including_meds_as_covariates\Transitions\Q3', sprintf('rt_trial_confirmatory2_meds_%s.csv', char(spec)));
        writetable(outT, out_csv);

        % Optional: make a forest figure from this CSV

        try
            out_fig = fullfile('D:\Papers\2025\Submitted\JAD (task microstates)\Analysis\Sensitivity tests\Including_meds_as_covariates\Transitions\Q3', sprintf('Fig_RT_Confirmatory2_meds_%s.tif', char(spec)));
            makeConfirmatory2Figure_fromFile(out_csv, out_fig);
        catch ME
            warning('Figure creation skipped for %s: %s', char(spec), ME.message);
        end

        fprintf('Q3 (meds=%s) saved: %s\n', char(spec), out_csv);
    end

end

if Link_trial_RT

    %% ================================================
    %  Link_trial_RT  — Reduced confirmatory family
    %  Trial-level RT ~ syntax with within-subject×condition demeaning
    %  Tests (Holm across 2):
    %    • N200 Anchor-2: stronger slope in valenced (Neg/Pos) vs Neutral?
    %    • P300 Anchor-2: stronger slope in valenced (Neg/Pos) vs Neutral?
    %  Negative control (not counted): LPP backbone
    %  Saves: rt_trial_confirmatory2.csv
    %  ================================================

    % -------------------------------------
    % Load demographics and RT
    % -------------------------------------

    load(fullfile(inputDir, 'France_data.mat'), 'France');

    % ---------------------------
    % Load outputStats(s)
    % ---------------------------

    S = load(fullfile(saveDir, 'outputStats.mat'));
    outputStats = S.outputStats;

    % -------------------------
    % ERP time windows
    % -------------------------

    windows.N200 = [180 300];
    windows.P300 = [300 500];
    windows.LPP = [500 1000];

    % ----------------------------------------------------------
    % HC topology → anchors and LPP backbone
    % ----------------------------------------------------------

    anchors = findAnchorsBackbone(outputStats, windows);

    % --------------------------------------------
    % Subject-ID map (order alignment)
    % --------------------------------------------

    subj_id_map = buildSubjectIdMap(France);

    % --------------------------
    % 1) Per-trial RT table
    % --------------------------

    rt_trials = buildRTTrialsTable(France); % subject, condition, trial, rt_ms, log_rt, age

    % -------------------------------------------------------
    % 2) Per-trial motif predictors from MSClass
    % -------------------------------------------------------

    trial_feats = computeTrialPredictors(outputStats, anchors, windows, subj_id_map);

    % ---------------------------------
    % 3) Join and prepare data
    % ---------------------------------

    trial_data = innerjoin(rt_trials, trial_feats, 'Keys', {'subject', 'condition', 'trial'});
    trial_data.condition = categorical(trial_data.condition, {'Neutral', 'Negative', 'Positive'});

    % Standardize predictors (z) over all trials

    pred_cols = {'N200_a1_rate', 'N200_a2_rate','P300_a1_rate', 'P300_a2_rate', 'LPP_backbone_rate'};

    for k = 1: numel(pred_cols)

        zname = ['z_' pred_cols{k}];
        mu = mean(trial_data.(pred_cols{k}), 'omitnan');
        sd = std(trial_data.(pred_cols{k}), 'omitnan');

        if sd == 0 || isnan(sd)
            trial_data.(zname) = zeros(height(trial_data), 1);
        else
            trial_data.(zname) = (trial_data.(pred_cols{k}) - mu) ./ sd;
        end
    end

    % Standardize trial index (time-on-task within condition)

    trial_data.trial_z = (trial_data.trial - mean(trial_data.trial)) ./ std(trial_data.trial);

    % ------------------------------------------
    % 4) Within subj×cond demeaning
    % ------------------------------------------

    pred_z = {'z_N200_a1_rate', 'z_N200_a2_rate', 'z_P300_a1_rate', 'z_P300_a2_rate', 'z_LPP_backbone_rate'};
    trial_data = demeanWithinSubjCond(trial_data, pred_z); % adds w_* (within) and m_* (cluster means)

    % Keep complete rows (only the three used predictors are strictly required)

    use = isfinite(trial_data.log_rt) & isfinite(double(trial_data.age)) & isfinite(trial_data.trial_z);
    need = {'z_N200_a2_rate', 'z_P300_a2_rate', 'z_LPP_backbone_rate'};

    for k = 1: numel(need)
        use = use & isfinite(trial_data.(['w_' need{k}])) & isfinite(trial_data.(['m_' need{k}]));
    end

    trialW = trial_data(use,:);

    % ----------------------------------------------------------
    % 5) Valence contrast (Neutral=0; Neg/Pos=1)
    % ----------------------------------------------------------

    trialW.c_val = double(trialW.condition ~= 'Neutral');

    % --------------------------------------------------------------
    % 6) Run the two confirmatory Anchor-2 models
    %    (random slope for within-effect; safe fallback)
    % --------------------------------------------------------------

    Rn200 = fit_one_anchor2(trialW, 'z_N200_a2_rate', 'N200 Anchor-2 (valenced vs neutral)');
    Rp300 = fit_one_anchor2(trialW, 'z_P300_a2_rate', 'P300 Anchor-2 (valenced vs neutral)');

    confirm2 = [Rn200; Rp300];

    % Holm step-down across the two p-values

    confirm2.p_Holm = holm_stepdown(confirm2.pValue);

    % -------------------------------------------------------
    % 7) Negative control (not counted in Holm)
    % -------------------------------------------------------

    Rctl = fit_one_anchor2(trialW, 'z_LPP_backbone_rate', 'LPP backbone (control; not counted)');
    Rctl.p_Holm = NaN;

    % -----------------------------
    % 8) Save + print results
    % -----------------------------

    outT = [confirm2; Rctl];

    if ~Pipeline_Sensitivity_tests
        writetable(outT, fullfile(saveDir, 'rt_trial_confirmatory2.csv'));
    else
        writetable(outT, fullfile(saveDir, 'rt_trial_confirmatory2_only_peaks.csv'));
    end

    disp('Reduced confirmatory family (interaction w*c_val = slope increase in valenced vs neutral):');
    disp(outT(:, {'Predictor', 'Beta', 'CI_low', 'CI_high', 'pValue', 'p_Holm', 'Beta_NeutralSlope', 'pValue_NeutralSlope', 'Beta_ValencedSlope', 'pValue_ValencedSlope'}));

    makeConfirmatory2Figure(saveDir);
end

function R = fit_one_anchor2(T, pred_base, pretty_name)

% Fit a single Anchor-2 (or control) model on within-demeaned data with valence contrast.
% Returns the interaction (confirmatory) and both simple slopes (Neutral and Valenced).

wcol = ['w_' pred_base]; % within-person deviation
mcol = ['m_' pred_base]; % Mundlak cluster mean (between-person)

% Keep complete rows for this predictor

use = isfinite(T.log_rt) & isfinite(double(T.age)) & isfinite(T.trial_z) & isfinite(T.(wcol)) & isfinite(T.(mcol)) & isfinite(T.c_val);
TT = T(use, :);

% Fixed effects and random effects (try slope, else fallback)

fix = ['log_rt ~ ' wcol ' + c_val + ' wcol ':c_val + ' mcol ' + age + trial_z'];
rf_full = ['(1 + ' wcol ' | subject)'];
rf_min = '(1 | subject)';

formula_try = [strrep(strrep(fix, 'wcol', wcol), 'mcol', mcol) ' + ' rf_full];
ok = true;

try
    lme = fitlme(TT, formula_try, 'FitMethod', 'REML');
catch
    ok = false;
end

if ~ok
    formula_try = [strrep(strrep(fix, 'wcol',wcol), 'mcol', mcol) ' + ' rf_min];
    lme = fitlme(TT, formula_try, 'FitMethod', 'REML');
end

C = lme.Coefficients; % table OR dataset
ncoef = size(C, 1); % works for both

% Names column exists in both forms

nam = string(C.Name);

% ---- Confirmatory coefficient = interaction (within × valence) ----

term_inter = string([wcol ':c_val']);
idxI = strcmp(nam, term_inter);

% Neutral slope (c_val = 0) = main effect of within

idxN = strcmp(nam, wcol);

% ------ Valenced slope = Neutral + Interaction (with CI & p) ------
% Point estimate:

b_val = C.Estimate(idxN) + C.Estimate(idxI);

% SE using covariance of fixed effects (guarded):

CovB = lme.CoefficientCovariance;
hasBoth = any(idxN) && any(idxI);

if hasBoth && size(CovB,1) >= ncoef && size(CovB,2) >= ncoef
    se_val = sqrt( CovB(idxN, idxN) + CovB(idxI, idxI) + 2 * CovB(idxN, idxI) );
else
    se_val = NaN;
end

% p-value for the linear combination (Neutral + Interaction)

L = zeros(1, ncoef);
L(idxN) = 1;
L(idxI) = 1;
p_val = coefTest(lme, L, 0); % returns p as first output

% 95% CI via t critical

vnames = getVarNamesLocal(C); % <— dataset/table safe
hasDFcol = any(strcmp(vnames, 'DF')) || any(strcmpi(vnames, 'DF'));

if hasDFcol && hasBoth

    % Access DF vector in a dataset-safe way

    if istable(C)
        dfN = C.DF(idxN);
        dfI = C.DF(idxI);
    else % dataset
        dfN = C.('DF')(idxN);
        dfI = C.('DF')(idxI);
    end

    if all(isfinite([dfN dfI]))
        df_eff = min([double(dfN) double(dfI)]);
    else
        df_eff = max(1, lme.DFE);
    end
else
    df_eff = max(1, lme.DFE);
end

if isfinite(se_val)
    tcrit = tinv(0.975, df_eff);
    ci_val = [b_val - tcrit * se_val, b_val + tcrit * se_val];
else
    ci_val = [NaN NaN];
end

% ---- Build output table ----

R = table;
R.Predictor = string(pretty_name);
R.Formula = string(formula_try);
R.Term = term_inter;

% Interaction (confirmatory)

R.Beta = C.Estimate(idxI);
R.CI_low = C.Lower(idxI);
R.CI_high = C.Upper(idxI);
R.pValue = C.pValue(idxI);

% Neutral slope

R.Beta_NeutralSlope = C.Estimate(idxN);
R.CI_low_NeutralSlope = C.Lower(idxN);
R.CI_high_NeutralSlope = C.Upper(idxN);
R.pValue_NeutralSlope = C.pValue(idxN);

% Valenced slope (Neutral + Interaction)

R.Beta_ValencedSlope = b_val;
R.CI_low_ValencedSlope = ci_val(1);
R.CI_high_ValencedSlope = ci_val(2);
R.pValue_ValencedSlope = p_val;

end

function q = holm_stepdown(p)

% Holm–Bonferroni step-down adjusted p-values for a vector p

p = p(:);
m = numel(p);
[ps, ix] = sort(p, 'ascend');
padj_sorted = zeros(m, 1);

for i = 1: m
    padj_sorted(i) = min(ps(i) * (m - i + 1), 1);

    if i > 1
        padj_sorted(i) = max(padj_sorted(i), padj_sorted(i - 1)); % monotone
    end
end

q = nan(size(p));
q(ix) = padj_sorted;

end

function anchors = findAnchorsBackbone(outputStats, windows)

% Identify HC-based anchors for N200/P300 and the strongest LPP reciprocal pair.

cond_list = {'negative', 'neutral', 'positive'};
counts_N200 = zeros(7, 7);
counts_P300 = zeros(7, 7);
counts_LPP = zeros(7, 7);

if ~isfield(outputStats, 'HC')
    error('HC field is missing in outputStats');
end

for c = 1: numel(cond_list)
    cond_name = cond_list{c};
    subj_arr = outputStats.HC.(cond_name);

    for s = 1: numel(subj_arr)
        X = subj_arr(s).MSClass

        if isempty(X)
            continue;
        end

        counts_N200 = counts_N200 + transitionCounts(X, windows.N200);
        counts_P300 = counts_P300 + transitionCounts(X, windows.P300);
        counts_LPP = counts_LPP + transitionCounts(X, windows.LPP);
    end
end

outdeg_N200 = sum(counts_N200, 2);
outdeg_P300 = sum(counts_P300, 2);
[~, idxN] = sort(outdeg_N200, 'descend');
[~, idxP] = sort(outdeg_P300, 'descend');
anchors.N200 = idxN(1: 2).';
anchors.P300 = idxP(1: 2).';
recip = counts_LPP + counts_LPP.';
recip(1: 8: 49) = 0;
maxVal = -inf;
pairIJ = [1 2];

for i = 1: 7
    for j = i + 1: 7
        val = recip(i, j);

        if val > maxVal
            maxVal = val;
            pairIJ = [i j];
        end
    end
end

anchors.LPP = pairIJ;

end

function C = transitionCounts(MSClass, win)

% 7×7 directed transition counts within samples win(1): win(2) across all trials

C = zeros(7, 7);
t1 = win(1);
t2 = win(2);
t1 = max(1, t1);
t2 = min(size(MSClass, 1), t2);

for tr = 1: size(MSClass, 2)
    l = double(MSClass(t1: t2, tr));
    valid = isfinite(l) & l >= 1 & l <= 7;
    l(~valid) = 0;
    pre = l(1: end - 1);
    nxt = l(2: end);
    mask = pre ~= nxt & pre > 0 & nxt > 0;
    pre = pre(mask);
    nxt = nxt(mask);

    for k = 1: numel(pre)
        C(pre(k), nxt(k)) = C(pre(k), nxt(k)) + 1;
    end
end
end

function subj_id_map = buildSubjectIdMap(France)

% Subject numbers in exact order of outputStats arrays

subj_id_map = struct;
subj_id_map.BP_I_Depressed = collect_ids_cell(France.BP_I.Depressed);
subj_id_map.BP_I_Euthymic = collect_ids_cell(France.BP_I.Euthymic);
subj_id_map.BP_II_Depressed = collect_ids_cell(France.BP_II.Depressed);
subj_id_map.BP_II_Euthymic = collect_ids_cell(France.BP_II.Euthymic);
subj_id_map.HC = collect_ids_cell(France.HC);
subj_id_map.Siblings = collect_ids_cell(France.Siblings);

    function ids = collect_ids_cell(cellarr)

        ids = NaN(1, numel(cellarr));

        for k = 1: numel(cellarr)
            s = cellarr{k};

            if isempty(s) || ~isfield(s, 'number')
                ids(k) = NaN;
            else
                ids(k) = s.number;
            end
        end
    end
end

function T = buildRTTrialsTable(France)

% Returns: subject, condition, trial, rt_ms, log_rt, age

rows = {};
add_group_mood(France, 'BP_I','Depressed');
add_group_mood(France, 'BP_I','Euthymic');
add_group_mood(France, 'BP_II','Depressed');
add_group_mood(France, 'BP_II','Euthymic');
add_group_nomood(France, 'HC');
add_group_nomood(France, 'Siblings');
T = vertcat(rows{:});

    function add_group_mood(F, gname, mname)

        if ~isfield(F, gname) || ~isfield(F.(gname), mname)
            return;
        end

        cells = F.(gname).(mname);

        for k = 1: numel(cells)
            s = cells{k};

            if isempty(s)
                continue;
            end

            add_three(s, 'Neutral');
            add_three(s, 'Negative');
            add_three(s, 'Positive');
        end
    end

    function add_group_nomood(F, gname)

        if ~isfield(F, gname)
            return;
        end

        cells = F.(gname);

        for k = 1: numel(cells)
            s = cells{k};

            if isempty(s)
                continue;
            end

            add_three(s, 'Neutral');
            add_three(s, 'Negative');
            add_three(s, 'Positive');
        end
    end

    function add_three(s, cond)

        v = s.(cond).Correct_resp_RT;
        v = v(:).';

        for t = 1: numel(v)
            rt = v(t);

            if ~isfinite(rt)
                continue;
            end

            R = table;
            R.subject = s.number;
            R.condition = string(cond);
            R.trial = t;
            R.rt_ms = rt;
            R.log_rt = log(rt);
            R.age = s.age;
            rows{end + 1, 1} = R;
        end
    end
end

function T = computeTrialPredictors(outputStats, anchors, windows, subj_id_map)

% Returns: subject, condition, trial, N200_a1_rate, N200_a2_rate, P300_a1_rate, P300_a2_rate, LPP_backbone_rate

rows = {};
conds = {'negative', 'neutral', 'positive'};
fields = fieldnames(outputStats);

for g = 1: numel(fields)
    grp = fields{g};

    if any(strcmpi(grp, {'HC', 'Siblings'}))
        ids_vec = subj_id_map.(grp);

        for c = 1: numel(conds)
            cname = conds{c};
            subj_arr = outputStats.(grp).(cname);

            for s = 1: numel(subj_arr)
                rec = subj_arr(s);
                subjnum = NaN;

                if s <= numel(ids_vec)
                    subjnum = ids_vec(s);
                end

                add_subject(rec, subjnum, cname);
            end
        end
    else

        ids_vec = subj_id_map.(grp);

        for c = 1: numel(conds)

            cname = conds{c};
            subj_arr = outputStats.(grp).(cname);

            for s = 1: numel(subj_arr)
                rec = subj_arr(s);
                subjnum = NaN;

                if s <= numel(ids_vec)
                    subjnum = ids_vec(s);
                end

                add_subject(rec, subjnum, cname);
            end

        end
    end
end

T = vertcat(rows{:});

    function add_subject(rec, subjnum, cname)

        X = rec.MSClass;

        if isempty(X)
            return;
        end

        nT = size(X, 2);
        L.N200 = (windows.N200(2) - windows.N200(1) + 1) / 1000;
        L.P300 = (windows.P300(2) - windows.P300(1) + 1) / 1000;
        L.LPP = (windows.LPP(2) - windows.LPP(1) +1) / 1000;

        for t = 1: nT

            Cn = transitionCounts_single(X(windows.N200(1): windows.N200(2), t));
            a1 = anchors.N200(1);
            a2 = anchors.N200(2);
            n200_a1 = (sum(Cn(a1, :)) - Cn(a1,a1)) / max(1e-9, L.N200);
            n200_a2 = (sum(Cn(a2, :)) - Cn(a2,a2)) / max(1e-9, L.N200);

            Cp = transitionCounts_single(X(windows.P300(1): windows.P300(2), t));
            p_a1 = anchors.P300(1);
            p_a2 = anchors.P300(2);
            p300_a1 = (sum(Cp(p_a1, :)) - Cp(p_a1, p_a1)) / max(1e-9, L.P300);
            p300_a2 = (sum(Cp(p_a2, :)) - Cp(p_a2, p_a2)) / max(1e-9, L.P300);

            Cl = transitionCounts_single(X(windows.LPP(1): windows.LPP(2), t));
            i = anchors.LPP(1);
            j = anchors.LPP(2);
            lpp_bb = (Cl(i, j) / max(1e-9, L.LPP) + Cl(j, i) / max(1e-9, L.LPP)) / 2;

            R = table;
            R.subject = subjnum;
            R.condition = string(upperFirst(cname));
            R.trial = t;
            R.N200_a1_rate = n200_a1;
            R.N200_a2_rate = n200_a2;
            R.P300_a1_rate = p300_a1;
            R.P300_a2_rate = p300_a2;
            R.LPP_backbone_rate = lpp_bb;
            rows{end + 1, 1} = R;
        end
    end
end

function C = transitionCounts_single(lbl)

% 7×7 transition counts for a single-trial label vector (time × 1)

C = zeros(7, 7);
l = double(lbl(:));
valid = isfinite(l) & l >= 1 & l <= 7;
l(~valid) = 0;
pre = l(1: end - 1);
nxt = l(2: end);
mask = pre ~= nxt & pre > 0 & nxt > 0;
pre = pre(mask);
nxt = nxt(mask);

for k = 1:numel(pre)
    C(pre(k), nxt(k)) = C(pre(k), nxt(k)) + 1;
end
end

function s = upperFirst(sin)

s = lower(string(sin));

if strlength(s) == 0
    return;
end

ch = char(s);
ch(1) = upper(ch(1));
s = string(ch);
end

function T = demeanWithinSubjCond(T, pred_cols)

% Create within-subject×condition centered (w_*) and cluster means (m_*)

[G, ~] = findgroups(T.subject, T.condition);

for k = 1: numel(pred_cols)
    col = pred_cols{k};
    mu_gc = splitapply(@nanmean, T.(col), G);
    T.(['m_' col]) = mu_gc(G);
    T.(['w_' col]) = T.(col) - T.(['m_' col]);
end
end

function vnames = getVarNamesLocal(C)

% Returns variable names for either table or dataset

if istable(C)
    vnames = C.Properties.VariableNames;
elseif isa(C, 'dataset')
    vnames = get(C, 'VarNames');
else
    vnames = {};
end
end

function outPath = makeConfirmatory2Figure(outputDir)

% Forest plot for reduced confirmatory family (Neutral vs Valenced slopes)
% Input:  outputDir containing rt_trial_confirmatory2.csv
% Output: saves Fig_RT_Confirmatory2.png in outputDir, returns its path

try
    csvPath = fullfile(outputDir, 'rt_trial_confirmatory2.csv');
catch
    csvPath = fullfile([outputDir '\Transition\Q3'], 'rt_trial_confirmatory2_only_peaks.csv');
end

if ~exist(csvPath, 'file')
    error('Could not find %s', csvPath);
end

T = readtable(csvPath);

% keep only the two Anchor-2 rows (drop LPP control if present)

keep = contains(T.Predictor, 'Anchor-2', 'IgnoreCase', true);
TT = T(keep, :);

% Order: P300 first, then N200 (top-to-bottom)

ord = {'P300', 'N200'};
rowOrder = zeros(height(TT), 1);

for i = 1: height(TT)

    if contains(TT.Predictor{i}, 'P300')
        rowOrder(i) = 1
    else
        rowOrder(i) = 2;
    end
end

[~, ix] = sort(rowOrder);
TT = TT(ix, :);

% Y positions

n = height(TT);
yBase = n: -1: 1;

% Pull quantities

bN = TT.Beta_NeutralSlope;
ciN = [TT.CI_low_NeutralSlope TT.CI_high_NeutralSlope];
bV = TT.Beta_ValencedSlope;
ciV = [TT.CI_low_ValencedSlope TT.CI_high_ValencedSlope];

% Aesthetics

colNeutral = [0.60 0.60 0.60];
colVal = [0.85 0.20 0.20];
mkN = 'o';
mkV = 's';
off = [-0.17, +0.17];

% Plot

fig = figure('Color', 'w', 'Position',[120 120 780 360]);
hold on;

for i = 1: n
    yN = yBase(i) + off(1);
    yV = yBase(i) + off(2);

    % Neutral CI + point

    plot(ciN(i,:), [yN yN], '-', 'Color', colNeutral, 'LineWidth', 2);
    plot(bN(i), yN, mkN, 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colNeutral);

    % Valenced CI + point (filled if Holm-significant interaction)

    plot(ciV(i, :), [yV yV], '-', 'Color', colVal, 'LineWidth', 2);
    fillVal = 'w';

    if ismember('p_Holm', TT.Properties.VariableNames) && isfinite(TT.p_Holm(i)) && TT.p_Holm(i) < 0.05
        fillVal = colVal;
    end

    plot(bV(i), yV, mkV, 'MarkerSize', 7, 'MarkerFaceColor', fillVal, 'MarkerEdgeColor', colVal);
end

xline(0, ':', 'Color', [0.4 0.4 0.4]);

% Y labels (clean names)

labs = strings(n, 1);

for i = 1: n
    if contains(TT.Predictor{i}, 'P300')
        labs(i) = "P300 Anchor-2";
    else
        labs(i) = "N200 Anchor-2";
    end
end

set(gca, 'YTick', sort(yBase), 'YTickLabel', labs, 'YDir','normal', 'Box', 'off');
xlabel('Standardized \beta on log RT (trial-level, within-subject×condition centered)');
legend({'Neutral','Valenced (Neg/Pos)'}, 'Location','southoutside','Orientation','horizontal');
title('Confirmatory behavior–syntax effects (Neutral vs Valenced slopes)\newline(Filled = Holm p<0.05 on interaction)');
grid on;

try
    outPath = fullfile(outputDir, 'Fig_RT_Confirmatory2.tif');
catch
    outPath = fullfile(outputDir, 'rt_trial_confirmatory2_only_peaks.tif');
end

exportgraphics(fig, outPath, 'Resolution', 300, 'BackgroundColor', 'white');
saveas(fig, outPath);
close(fig);

fprintf('Saved %s\n', outPath);
end

function medsT = build_meds_subject_table(France, spec)

% Returns a table with one row per subject: subject, [AD AP MS ANX OTHER] or [poly_count]
% spec: 'classes' | 'poly' (defaults to 'classes' if empty)

if nargin < 2 || isempty(spec)
    spec = 'classes';
end

spec = validatestring(spec, {'classes', 'poly'});

rows = {};
add_group_mood(France, 'BP_I', 'Depressed');
add_group_mood(France, 'BP_I', 'Euthymic');
add_group_mood(France, 'BP_II', 'Depressed');
add_group_mood(France, 'BP_II', 'Euthymic');
add_group_nomood(France, 'HC');
add_group_nomood(France, 'Siblings');

medsT = vertcat(rows{:});

% Deduplicate (just in case) and keep one row per subject

[~, ix] = unique(medsT.subject);
medsT = medsT(sort(ix), :);

    function add_group_mood(F, gname, mname)

        if ~isfield(F,gname) || ~isfield(F.(gname), mname)
            return;
        end

        cells = F.(gname).(mname);

        for k = 1: numel(cells)

            s = cells{k};
            if isempty(s) || ~isfield(s, 'number')
                continue;
            end

            rows{end + 1, 1} = one_row(s);
        end
    end

    function add_group_nomood(F,gname)

        if ~isfield(F, gname)
            return;
        end

        cells = F.(gname);

        for k = 1: numel(cells)
            s = cells{k};

            if isempty(s) || ~isfield(s, 'number')
                continue;
            end

            rows{end + 1, 1} = one_row(s);
        end
    end

    function T = one_row(s)

        subject = double(s.number);

        % Safe pulls with defaults (missing -> 0)

        if ~isfield(s,'med') || isempty(s.med)
            M = struct();
        else
            M = s.med;
        end
        getf = @(f,def) double((isfield(M,f) && ~isempty(M.(f))) * M.(f) + (~isfield(M, f) || isempty(M.(f))) * def);

        switch spec
            case 'classes'
                AD = getf('AD', 0);
                AP = getf('AP', 0);
                MS = getf('MS', 0);
                ANX = getf('ANX', 0);
                OTHER = getf('OTHER', 0);
                T = table(subject, AD, AP, MS, ANX, OTHER);

            case 'poly'
                poly_count = getf('poly_count', 0);
                T = table(subject, poly_count);
        end
    end
end

function R = fit_one_anchor2_meds(T, pred_base, pretty_name, medVars)

% As fit_one_anchor2, but includes medication covariates (main effects).
% Inputs:
% T  table with trialW + meds columns
%  pred_base  e.g., 'z_N200_a2_rate'
%  pretty_name label for reporting
%  medVars  cellstr/string array of med column names to include

if nargin < 4
    medVars = strings(0, 1);
end

medVars = string(medVars(:))';

wcol = ['w_' pred_base]; % within-person deviation
mcol = ['m_' pred_base]; % Mundlak cluster mean

% Keep complete rows for this predictor + meds

use = isfinite(T.log_rt) & isfinite(double(T.age)) & isfinite(T.trial_z) & isfinite(T.(wcol)) & isfinite(T.(mcol)) & isfinite(T.c_val);

for v = medVars
    use = use & isfinite(T.(v));
end

TT = T(use, :);

% Fixed effects (meds enter additively as covariates)

if isempty(medVars)
    medTerm = '';
else
    medTerm = [' + ' strjoin(cellstr(medVars), ' + ')];
end

fix = ['log_rt ~ ' wcol ' + c_val + ' wcol ':c_val + ' mcol ' + age + trial_z' medTerm];
rf_full = ['(1 + ' wcol ' | subject)'];
rf_min = '(1 | subject)';

formula_try = [fix ' + ' rf_full];
ok = true;

try
    lme = fitlme(TT, formula_try, 'FitMethod', 'REML');
catch
    ok = false;
end

if ~ok
    formula_try = [fix ' + ' rf_min];
    lme = fitlme(TT, formula_try, 'FitMethod', 'REML');
end

C = lme.Coefficients;
nam = string(C.Name);
ncoef = size(C, 1);

% Confirmatory coefficient = interaction (within × valence)

term_inter = string([wcol ':c_val']);
idxI = strcmp(nam, term_inter);

% Neutral slope = main effect of within

idxN = strcmp(nam, wcol);

% Valenced slope = Neutral + Interaction (with CI & p via linear combo)

b_val = C.Estimate(idxN) + C.Estimate(idxI);
CovB = lme.CoefficientCovariance;
hasBoth = any(idxN) && any(idxI);

if hasBoth && size(CovB, 1) >= ncoef && size(CovB, 2) >= ncoef
    se_val = sqrt( CovB(idxN, idxN) + CovB(idxI, idxI) + 2 * CovB(idxN, idxI) );
else
    se_val = NaN;
end

L = zeros(1, ncoef);
L(idxN) = 1;
L(idxI) = 1;
p_val = coefTest(lme, L, 0); % p for Valenced simple slope

% DF for CI

vnames = getVarNamesLocal(C);
hasDFcol = any(strcmp(vnames, 'DF')) || any(strcmpi(vnames, 'DF'));

if hasDFcol && hasBoth
    if istable(C)
        dfN = C.DF(idxN);
        dfI = C.DF(idxI);
    else
        dfN = C.('DF')(idxN);
        dfI = C.('DF')(idxI);
    end

    if all(isfinite([dfN dfI]))
        df_eff = min([double(dfN) double(dfI)]);
    else
        df_eff = max(1, lme.DFE);
    end
else
    df_eff = max(1, lme.DFE);
end

if isfinite(se_val)
    tcrit = tinv(0.975, df_eff);
    ci_val = [b_val - tcrit * se_val, b_val + tcrit * se_val];
else
    ci_val = [NaN NaN];
end

% Assemble output (same columns as baseline helper)

R = table;
R.Predictor = string(pretty_name);
R.Formula = string(formula_try);
R.Term = term_inter;

R.Beta = C.Estimate(idxI);
R.CI_low = C.Lower(idxI);
R.CI_high = C.Upper(idxI);
R.pValue = C.pValue(idxI);

R.Beta_NeutralSlope = C.Estimate(idxN);
R.CI_low_NeutralSlope = C.Lower(idxN);
R.CI_high_NeutralSlope = C.Upper(idxN);
R.pValue_NeutralSlope = C.pValue(idxN);

R.Beta_ValencedSlope = b_val;
R.CI_low_ValencedSlope = ci_val(1);
R.CI_high_ValencedSlope = ci_val(2);
R.pValue_ValencedSlope = p_val;
end

function makeConfirmatory2Figure_fromFile(csvPath, outPath)

% Forest plot for the two Anchor-2 models (reads CSV you just wrote)

T = readtable(csvPath);
keep = contains(T.Predictor, 'Anchor-2', 'IgnoreCase', true);
TT = T(keep, :);

% Order: P300 first, then N200

rowOrder = zeros(height(TT), 1);

for i = 1: height(TT)
    if contains(TT.Predictor{i}, 'P300')
        rowOrder(i) = 1;
    else
        rowOrder(i) = 2;
    end
end

[~, ix] = sort(rowOrder);
TT = TT(ix, :);
n = height(TT);
yBase = n: -1: 1;

bN = TT.Beta_NeutralSlope;
ciN = [TT.CI_low_NeutralSlope TT.CI_high_NeutralSlope];
bV = TT.Beta_ValencedSlope;
ciV = [TT.CI_low_ValencedSlope TT.CI_high_ValencedSlope];

colNeutral = [0.60 0.60 0.60];
colVal = [0.85 0.20 0.20];
mkN = 'o';
mkV = 's';
off = [-0.17, +0.17];

fig = figure('Color','w','Position',[120 120 780 360]);
hold on;

for i = 1: n
    yN = yBase(i) + off(1);
    yV = yBase(i) + off(2);
    plot(ciN(i, :), [yN yN], '-', 'Color', colNeutral, 'LineWidth', 2);
    plot(bN(i), yN, mkN, 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colNeutral);

    plot(ciV(i, :), [yV yV], '-', 'Color', colVal, 'LineWidth', 2);
    fillVal = 'w';

    if ismember('p_Holm', TT.Properties.VariableNames) && isfinite(TT.p_Holm(i)) && TT.p_Holm(i) < 0.05
        fillVal = colVal;
    end
    plot(bV(i), yV, mkV, 'MarkerSize', 7, 'MarkerFaceColor', fillVal, 'MarkerEdgeColor', colVal);
end

xline(0, ':', 'Color', [0.4 0.4 0.4]);

labs = strings(n, 1);

for i = 1: n
    if contains(TT.Predictor{i}, 'P300')
        labs(i) = "P300 Anchor-2";
    else
        labs(i) = "N200 Anchor-2";
    end
end

set(gca, 'YTick', sort(yBase), 'YTickLabel', labs, 'YDir', 'normal', 'Box', 'off');
xlabel('Standardized \beta on log RT (trial-level, within-subject×condition centered)');
legend({'Neutral', 'Valenced (Neg/Pos)'}, 'Location','southoutside', 'Orientation',' horizontal');
title('Confirmatory behavior–syntax effects with medication controls\newline(Filled = Holm p<0.05 on interaction)');
grid on;

if nargin < 2 || isempty(outPath)
    [folder,~,~] = fileparts(csvPath);
    outPath = fullfile(folder, 'Fig_RT_Confirmatory2_meds.tif');
end

exportgraphics(fig, outPath, 'Resolution', 300, 'BackgroundColor', 'white');
close(fig);
fprintf('Saved %s\n', outPath);

end