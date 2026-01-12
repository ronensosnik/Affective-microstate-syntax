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

%%

% ==== First stage: Extract the microstates. Next, compute the three static measures ========

Find_microstates = 0;
Compute_and_save_metrics = 0;

% ==== Second stage: Analyze static metrics and transitions (Q1 and Q2)  ========

Analyze_metrics = 0;
Analyze_transitions = 1;

% ==== Third stage: test for a link between transitions and RT  ========

Link_trial_RT = 0;

% ===== Fourth stage: Sensitivity tests: a) Smoothness penalty, b) Fit peaks only ======

Peaks_only_sensitivity_test = 0;
Smootness_penalty_sensitivity_test = 0;

% =========================================================================================

inputDir = 'D:\Papers\2025\In preparation\XXX (task microstates)\Data\';

if ~Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
    saveDir = 'D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\';
    ClustPar_GFPPeaks = false; FitPar_PeakFit = 0;
elseif Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
    saveDir = 'D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\Sensitivity tests\GFP_Peaks only\';
    ClustPar_GFPPeaks = true; FitPar_PeakFit = 1;
elseif ~Peaks_only_sensitivity_test && Smootness_penalty_sensitivity_test
    saveDir = 'D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\Sensitivity tests\Smoothness_penalty\';
    ClustPar_GFPPeaks = false; FitPar_PeakFit = 0;
end

ERPs = {'N200', 'P300', 'LPP'};
Microstates = {'Microstate_A', 'Microstate_B', 'Microstate_C', 'Microstate_D', 'Microstate_E', 'Microstate_F', 'Microstate_G'};
Groups = {'BP_I_Depressed', 'BP_I_Euthymic', 'BP_II_Depressed', 'BP_II_Euthymic', 'HC', 'Siblings'};
Conditions = {'Negative', 'Neutral', 'Positive'};
Occurrence = {'Occurrence_A', 'Occurrence_B', 'Occurrence_C' ,'Occurrence_D' ,'Occurrence_E' ,'Occurrence_F', 'Occurrence_G'};

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
    Backfit_and_quantify_temp_dynam = 0;

    %% Setting all parameters

    groupFolders = dir(inputDir);
    groupFolders = groupFolders([groupFolders.isdir]);
    groupFolders = groupFolders(~matches({groupFolders.name}, {'.', '..'}));
    nGroups = length(groupFolders);
    groupNames = cell(1, nGroups);
    dataDirs = cell(1, nGroups);

    % Set clustering parameters

    ClustPar.UseAAHC = false;
    ClustPar.MinClasses = 7;
    ClustPar.MaxClasses = 7;
    ClustPar.MaxMaps = inf;
    ClustPar.GFPPeaks = ClustPar_GFPPeaks;
    ClustPar.IgnorePolarity = true;
    ClustPar.Normalize = true;
    ClustPar.Allow_early_stop = true;

    ClustPar.Restarts = 100;

    % Set backfitting parameters

    FitPar.Classes = 7;
    FitPar.PeakFit = FitPar_PeakFit;

    if ~Smootness_penalty_sensitivity_test
        FitPar.lambda = 0.3;
    else
        FitPar.lambda = 0.7;
    end

    FitPar.b = 20;

    for i = 1: nGroups
        groupDir = fullfile(inputDir, groupFolders(i).name);
        groupNames{i} = groupFolders(i).name;
        subjFolders = dir(groupDir);
        subjFolders = subjFolders(~matches({subjFolders.name}, {'.', ' ..'}));
        subjNames{i} = {subjFolders.name};
        dataDirs{i} = cellfun(@(x) fullfile(groupDir, x), subjNames{i}, 'UniformOutput', false);
    end

    % Start EEGLAB and find Microstates plugin files

    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
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
        save([saveDir ['7_clusters_' ERP_to_use_for_clustering.name '.mat']], 'ALLEEG', 'EEG', 'CURRENTSET', 'Number_of_subjects', 'GroupIdx', '-v7.3');

        loaded = 1;

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

        if ~loaded
            load([saveDir ['\7_clusters_' ERP_to_use_for_clustering.name '.mat']]);
            loaded = 1;
        end

        pop_ShowIndMSMaps(ALLEEG, Number_of_subjects + 1: Number_of_subjects + 7, 'Classes', 7);

        load('D:\eeglab2025.1.0\plugins\MICROSTATELAB2.1\Templates\Custo.mat');
        pop_ShowIndMSMaps(Custo, 1, 'Classes', 7);
    end

    %% Compare classes across groups

    if Compare_classes_across_groups

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

        if ~loaded
            load([saveDir ['\7_clusters_' ERP_to_use_for_clustering.name '.mat']]);
            loaded = 1;
        end

        % Compare grand mean and Custo et al template

        sharedVarTable_GrandMean = pop_CompareMSMaps(ALLEEG, 'MeanSets', 79, 'PublishedSets', 'Custo2017', 'Classes', 7, 'Filename', [saveDir '\sharedvars.mat'], 'gui', 0);
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Within each window, labels were clipped to the window bounds.
    % For each microstate class we compute:
    % (i) coverage, the proportion of samples in the window carrying that label (segments spanning a boundary contributed by their overlap only);
    % (ii) occurrence, the number of onsets whose first sample fell within the window, divided by the window length to yield a rate (Hz);
    % (iii) mean duration, the average duration of segments whose onsets occurred within the window (segments truncated by the window boundary were excluded from the duration calculation).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    load('D:\Papers\2025\In preparation\XXX (task microstates)\Data\France_data.mat');
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

                for t = 1: numel(Conditions)
                    c = Conditions{t};
                    Data.(e).(m).(g).(c) = repmat({leaf}, 1, 1); % create the 6 fields here
                    Mean.(e).(m).(g).(c) = repmat({leaf}, 1, 1); % create the 6 fields here
                end
            end
        end
    end

    for erp = 1: 3
        for micro = 1: 7
            for gro = 1: 6
                for condd = 1: 3

                    clear vvv; vvv = {outputStats.(Groups{gro}).([lower(Conditions{condd}(1)) Conditions{condd}(2: end)]).MSClass};

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
                                Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.num_occurrences = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.num_occurrences aaa(micro) * (1000 / length(ERP.window{erp}))];
                            catch
                                Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.num_occurrences = [];
                                Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.num_occurrences = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.num_occurrences aaa(micro) * (1000 / length(ERP.window{erp}))];
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
                                Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.duration = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.duration mean(z_ee(y_ee == micro))];

                            catch
                                Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.duration = [];
                                Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.duration = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.duration mean(z_ee(y_ee == micro))];
                            end

                            % coverage

                            try
                                Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.coverage = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.coverage (sum(z(y == micro))) / (ERP.window{erp}(end) - ERP.window{erp}(1) + 1)];
                            catch
                                Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.coverage = [];
                                Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.coverage = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.coverage (sum(z(y == micro))) / (ERP.window{erp}(end) - ERP.window{erp}(1) + 1)];
                            end

                            % subject number, age and medications

                            switch Groups{gro}
                                case 'BP_I_Depressed'
                                    try
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age France.BP_I.Depressed{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.BP_I.Depressed{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds France.BP_I.Depressed{subj_num}.med];
                                    catch
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age France.BP_I.Depressed{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.BP_I.Depressed{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds France.BP_I.Depressed{subj_num}.med];
                                    end

                                case 'BP_II_Depressed'
                                    try
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age France.BP_II.Depressed{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.BP_II.Depressed{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.BP_II.Depressed{subj_num}.med];
                                    catch
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age France.BP_II.Depressed{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.BP_II.Depressed{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds France.BP_II.Depressed{subj_num}.med];
                                    end

                                case 'BP_I_Euthymic'
                                    try
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age France.BP_I.Euthymic{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.BP_I.Euthymic{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds France.BP_I.Euthymic{subj_num}.med];
                                    catch
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age France.BP_I.Euthymic{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.BP_I.Euthymic{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds France.BP_I.Euthymic{subj_num}.med];
                                    end

                                case 'BP_II_Euthymic'
                                    try
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age France.BP_II.Euthymic{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.BP_II.Euthymic{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.BP_II.Euthymic{subj_num}.med];
                                    catch
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age France.BP_II.Euthymic{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.BP_II.Euthymic{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds France.BP_II.Euthymic{subj_num}.med];
                                    end

                                case 'HC'
                                    try
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age France.HC{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.HC{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.HC{subj_num}.med];
                                    catch
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age France.HC{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.HC{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds France.HC{subj_num}.med];
                                    end

                                case 'Siblings'
                                    try
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age France.Siblings{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.Siblings{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.Siblings{subj_num}.med];
                                    catch
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age France.Siblings{subj_num}.age];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number France.Siblings{subj_num}.number];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [];
                                        Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = [Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds France.Siblings{subj_num}.med];
                                    end

                            end
                        end

                    end

                    for subj_num = 1: size(vvv, 2)

                        Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.num_occurrences = mean(Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.num_occurrences);
                        Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.duration = nanmean(Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.duration);
                        Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.coverage = mean(Data.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.coverage) * 100;
                        Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.group = Groups{gro};
                        Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.condition = Conditions{condd};

                        switch Groups{gro}

                            case 'BP_I_Depressed'
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = France.BP_I.Depressed{subj_num}.number;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = France.BP_I.Depressed{subj_num}.age;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = France.BP_I.Depressed{subj_num}.med;
                            case 'BP_II_Depressed'
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = France.BP_II.Depressed{subj_num}.number;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = France.BP_II.Depressed{subj_num}.age;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = France.BP_II.Depressed{subj_num}.med;
                            case 'BP_I_Euthymic'
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = France.BP_I.Euthymic{subj_num}.number;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = France.BP_I.Euthymic{subj_num}.age;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = France.BP_I.Euthymic{subj_num}.med;
                            case 'BP_II_Euthymic'
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = France.BP_II.Euthymic{subj_num}.number;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = France.BP_II.Euthymic{subj_num}.age;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = France.BP_II.Euthymic{subj_num}.med;
                            case 'HC'
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = France.HC{subj_num}.number;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = France.HC{subj_num}.age;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = France.HC{subj_num}.med;
                            case 'Siblings'
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.number = France.Siblings{subj_num}.number;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.age = France.Siblings{subj_num}.age;
                                Mean.(ERPs{erp}).(Microstates{micro}).(Groups{gro}).(Conditions{condd}){subj_num}.meds = France.Siblings{subj_num}.med;
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

            ALL_Mean.(eName).(mName).Age = [];
            ALL_Mean.(eName).(mName).Number = [];
            ALL_Mean.(eName).(mName).Meds = [];
            ALL_Mean.(eName).(mName).Coverage = [];
            ALL_Mean.(eName).(mName).Duration = [];
            ALL_Mean.(eName).(mName).Num_occurrence = [];
            ALL_Mean.(eName).(mName).Group = {};
            ALL_Mean.(eName).(mName).Condition = {};

            for gro = 1: 6

                gName = Groups{gro};

                for condd = 1: 3

                    cName = Conditions{condd};

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
                    Conds = cellfun(@(s) s.condition, S, 'UniformOutput', false);

                    % append (concatenate) to accumulators

                    ALL_Mean.(eName).(mName).Age = [ALL_Mean.(eName).(mName).Age, ages];
                    ALL_Mean.(eName).(mName).Meds = [ALL_Mean.(eName).(mName).Meds, meds];
                    ALL_Mean.(eName).(mName).Number = [ALL_Mean.(eName).(mName).Number, nums];
                    ALL_Mean.(eName).(mName).Coverage = [ALL_Mean.(eName).(mName).Coverage, covs];
                    ALL_Mean.(eName).(mName).Duration = [ALL_Mean.(eName).(mName).Duration, durs];
                    ALL_Mean.(eName).(mName).Num_occurrence = [ALL_Mean.(eName).(mName).Num_occurrence, occs];
                    ALL_Mean.(eName).(mName).Group = [ALL_Mean.(eName).(mName).Group, grps];
                    ALL_Mean.(eName).(mName).Condition = [ALL_Mean.(eName).(mName).Condition, Conds];
                end
            end
        end
    end

    save(fullfile(saveDir, ['Metrics_ALL.mat']), 'ALL_Mean');

end

if Analyze_metrics

    load([saveDir ['\Metrics_ALL.mat']]);

    % For each of the three measures, we run 21 sets of models (7 microstates * 3 ERPs). We apply ERP BH-FDR correction.

    if ~Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
        Out_Duration_base = lme_posthoc(ALL_Mean, 'Duration');
        Out_Occurrence_base = lme_posthoc(ALL_Mean, 'Num_occurrence');
        Out_Coverage_base = lme_posthoc(ALL_Mean, 'Coverage');
    elseif Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
        Out_Duration_sensitivity_peaks = lme_posthoc(ALL_Mean, 'Duration');
        Out_Occurrence_sensitivity_peaks = lme_posthoc(ALL_Mean, 'Num_occurrence');
        Out_Coverage_sensitivity_peaks = lme_posthoc(ALL_Mean, 'Coverage');
    elseif ~Peaks_only_sensitivity_test && Smootness_penalty_sensitivity_test
        Out_Duration_sensitivity_smoothness = lme_posthoc(ALL_Mean, 'Duration');
        Out_Occurrence_sensitivity_smoothness = lme_posthoc(ALL_Mean, 'Num_occurrence');
        Out_Coverage_sensitivity_smoothness = lme_posthoc(ALL_Mean, 'Coverage');
    end

    % Report

    if ~Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
        print_summary_measure(Out_Duration_base);
        print_summary_measure(Out_Occurrence_base);
        print_summary_measure(Out_Coverage_base);
    elseif Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
        print_summary_measure(Out_Duration_sensitivity_peaks);
        print_summary_measure(Out_Occurrence_sensitivity_peaks);
        print_summary_measure(Out_Coverage_sensitivity_peaks);
    elseif ~Peaks_only_sensitivity_test && Smootness_penalty_sensitivity_test
        print_summary_measure(Out_Duration_sensitivity_smoothness);
        print_summary_measure(Out_Occurrence_sensitivity_smoothness);
        print_summary_measure(Out_Coverage_sensitivity_smoothness);
    end

    % Generate plots for any ERP BH-FDR-significant effects

    if ~Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
        plot_emm_significant(Out_Duration_base, 'Duration', ALL_Mean);
        plot_emm_significant(Out_Occurrence_base, 'Num_occurrence', ALL_Mean);
        plot_emm_significant(Out_Coverage_base, 'Coverage', ALL_Mean);
    end

    % Save results in .mat file and export to .csv

    if ~Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test

        save(fullfile(saveDir, ['Out_Duration_base.mat']), 'Out_Duration_base');
        save(fullfile(saveDir, ['Out_Occurrence_base.mat']), 'Out_Occurrence_base');
        save(fullfile(saveDir, ['Out_Coverage_base.mat']), 'Out_Coverage_base');

        Export_pairwise(Out_Duration_base, 'Duration', 'baseline', fullfile(saveDir, ['Out_Duration_baseline.csv']));
        Export_pairwise(Out_Occurrence_base, 'Occurrence', 'baseline', fullfile(saveDir, ['Out_Occurrence_base.csv']));
        Export_pairwise(Out_Coverage_base, 'Coverage', 'baseline', fullfile(saveDir, ['Out_Coverage_base.csv']));

    elseif Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test

        save(fullfile([saveDir '\Static'], ['Out_Duration_sensitivity_peaks.mat']), 'Out_Duration_sensitivity_peaks');
        save(fullfile([saveDir '\Static'], ['Out_Occurrence_sensitivity_peaks.mat']), 'Out_Occurrence_sensitivity_peaks');
        save(fullfile([saveDir '\Static'], ['Out_Coverage_sensitivity_peaks.mat']), 'Out_Coverage_sensitivity_peaks');

        Export_pairwise(Out_Duration_sensitivity_peaks, 'Duration', 'baseline', fullfile([saveDir '\Static'], ['Out_Duration_sensitivity_peaks.csv']));
        Export_pairwise(Out_Occurrence_sensitivity_peaks, 'Occurrence', 'baseline', fullfile([saveDir '\Static'], ['Out_Occurrence_sensitivity_peaks.csv']));
        Export_pairwise(Out_Coverage_sensitivity_peaks, 'Coverage', 'baseline', fullfile([saveDir '\Static'], ['Out_Coverage_sensitivity_peaks.csv']));

        % Compare baseline pipeline with smoothness panalty pipeline

        load('D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\Out_Duration_base.mat', 'Out_Duration_base');
        load('D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\Out_Occurrence_base.mat', 'Out_Occurrence_base');
        load('D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\Out_Coverage_base.mat', 'Out_Coverage_base');

        Duration_baseline_vs_sensitivity_peaks = summarize_metric_sensitivity_robustness(Out_Duration_base, Out_Duration_sensitivity_peaks, 'Duration');
        Occurrence_baseline_vs_sensitivity_peaks = summarize_metric_sensitivity_robustness(Out_Occurrence_base, Out_Occurrence_sensitivity_peaks, 'Occurrence');
        Coverage_baseline_vs_sensitivity_peaks = summarize_metric_sensitivity_robustness(Out_Coverage_base, Out_Coverage_sensitivity_peaks, 'Coverage');

        save(fullfile([saveDir '\Static'], ['Duration_baseline_vs_sensitivity_peaks.mat']), 'Duration_baseline_vs_sensitivity_peaks');
        save(fullfile([saveDir '\Static'], ['Occurrence_baseline_vs_sensitivity_peaks.mat']), 'Occurrence_baseline_vs_sensitivity_peaks');
        save(fullfile([saveDir '\Static'], ['Coverage_baseline_vs_sensitivity_peaks.mat']), 'Coverage_baseline_vs_sensitivity_peaks');

    elseif ~Peaks_only_sensitivity_test && Smootness_penalty_sensitivity_test

        save(fullfile([saveDir '\Static'], ['Out_Duration_sensitivity_smoothness.mat']), 'Out_Duration_sensitivity_smoothness');
        save(fullfile([saveDir '\Static'], ['Out_Occurrence_sensitivity_smoothness.mat']), 'Out_Occurrence_sensitivity_smoothness');
        save(fullfile([saveDir '\Static'], ['Out_Coverage_sensitivity_smoothness.mat']), 'Out_Coverage_sensitivity_smoothness');

        Export_pairwise(Out_Duration_sensitivity_smoothness, 'Duration', 'baseline', fullfile([saveDir '\Static'], ['Out_Duration_sensitivity_smoothness.csv']));
        Export_pairwise(Out_Occurrence_sensitivity_smoothness, 'Occurrence', 'baseline', fullfile([saveDir '\Static'], ['Out_Occurrence_sensitivity_smoothness.csv']));
        Export_pairwise(Out_Coverage_sensitivity_smoothness, 'Coverage', 'baseline', fullfile([saveDir '\Static'], ['Out_Coverage_sensitivity_smoothness.csv']));

        % Compare baseline pipeline with smoothness panalty pipeline

        load('D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\Out_Duration_base.mat', 'Out_Duration_base');
        load('D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\Out_Occurrence_base.mat', 'Out_Occurrence_base');
        load('D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\Out_Coverage_base.mat', 'Out_Coverage_base');

        Duration_baseline_vs_sensitivity_robustness = summarize_metric_sensitivity_robustness(Out_Duration_base, Out_Duration_sensitivity_smoothness, 'Duration');
        Occurrence_baseline_vs_sensitivity_robustness = summarize_metric_sensitivity_robustness(Out_Occurrence_base, Out_Occurrence_sensitivity_smoothness, 'Occurrence');
        Coverage_baseline_vs_sensitivity_robustness = summarize_metric_sensitivity_robustness(Out_Coverage_base, Out_Coverage_sensitivity_smoothness, 'Coverage');

        save(fullfile([saveDir '\Static'], ['Duration_baseline_vs_sensitivity_robustness.mat']), 'Duration_baseline_vs_sensitivity_robustness');
        save(fullfile([saveDir '\Static'], ['Occurrence_baseline_vs_sensitivity_robustness.mat']), 'Occurrence_baseline_vs_sensitivity_robustness');
        save(fullfile([saveDir '\Static'], ['Coverage_baseline_vs_sensitivity_robustness.mat']), 'Coverage_baseline_vs_sensitivity_robustness');
    end

end

if Analyze_transitions

    %% ========================================================================= %%
    %  TASK MICROSTATES — TWO TRANSITION ANALYSES (AGE-ADJUSTED, BH-FDR PER ERP)
    %  Q1. Within each panel (ERP × Group × Condition), which X -> Y transitions
    %         deviate from independence? [Poisson GLMs with independence offset]
    %         Multiplicity: BH-FDR pooled per ERP across all Groups × Conditions × 42 edges.
    %
    %  Q2. Anchored contrasts (Target vs HC at same age) handled in buildAnchoredContrasts.m
    %         Multiplicity: BH-FDR pooled per ERP across all Targets × Conditions × 42 edges.
    %% ========================================================================== %%

    %% ---------------------------------------- I/O ------------------------------------------- %%

    data_path = 'D:\Papers\2025\In preparation\XXX (task microstates)\Data\France_data.mat';
    stats_path = [saveDir '\outputStats.mat'];

    load(data_path, 'France');
    load(stats_path, 'outputStats');

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
    % trans_counts (7 x 7; NaN diag), M_segments, N_transitions, C_segments (1 x 7), Occurrence_*.

    Out = struct();
    MicrostateLabels = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};

    for e = 1: numel(ERPs)
        erp = ERPs{e};

        for c = 1: numel(Conditions)
            condField = lower(Conditions{c});

            for g = 1: numel(Groups)
                grp = Groups{g};
                subj_list = outputStats.(grp).(condField);
                nSubj = numel(subj_list);

                if isfield(AgeMap, grp) && numel(AgeMap.(grp)) ~= nSubj
                    warning('%s: Age vector length (%d) != #subjects (%d).', grp, numel(AgeMap.(grp)), nSubj);
                end


                for sj = 1: nSubj
                    S = subj_list(sj);

                    if ~isfield(S, 'MSClass') || isempty(S.MSClass)
                        continue;
                    end

                    addTrial = ones(7, 7);
                    Trans_dur = repmat({[nan nan]}, 7, 7);

                    segments = cell(1, size(S.MSClass, 2));
                    Duration = cell(1, size(S.MSClass, 2));
                    num_transitions = zeros(1, size(S.MSClass, 2));
                    Transition_matrix = zeros(7, 7);
                    Transition_matrix(1: 8: 49) = NaN;

                    for tri = 1: size(S.MSClass, 2)

                        aa = [S.MSClass(ERP.window{e}, tri); -999];
                        clear nonrepeats; [nonrepeats, Dura] = find(diff(aa) ~= 0);
                        segments{tri} = aa(nonrepeats).';
                        Duration{tri}(1) = nonrepeats(1);

                        try
                            Duration{tri}(2: length(nonrepeats)) = diff(nonrepeats)'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        catch
                        end

                        clear aaaa bbbb; [aaaa, bbbb] = find(segments{tri} == 0);

                        if ~isempty(bbbb)
                            Duration{tri}(bbbb) = 0;
                        end

                        segments{tri} = segments{tri}(segments{tri} ~= 0);
                        Duration{tri} = Duration{tri}(Duration{tri} ~= 0);

                        num_transitions(tri) = max(0, numel(segments{tri}) - 1);

                        for p = 1: (numel(segments{tri}) - 1)

                            i = segments{tri}(p);
                            j = segments{tri}(p + 1);

                            if i ~= j

                                Transition_matrix(i, j) = Transition_matrix(i, j) + 1;

                                if p > 1 && (p + 1) < numel(segments{tri}) % Collect the duration of the source and target microstates only if they started and ended within the trial

                                    if ~isnan(Trans_dur{i, j}(1))
                                        addTrial(i, j) = addTrial(i, j) + 1;
                                    end

                                    Trans_dur{i, j}(addTrial(i, j), 1) = Duration{tri}(p);
                                    Trans_dur{i, j}(addTrial(i, j), 2) = Duration{tri}(p + 1);
                                end
                            end
                        end
                    end

                    clear M_segments; M_segments = sum(cellfun(@numel, segments));
                    clear pool; pool = [segments{:}];
                    C_segments = zeros(1, 7);

                    for s = 1: 7
                        C_segments(s) = nnz(pool == s);
                    end

                    clear N_transitions; N_transitions = sum(num_transitions);

                    Out.(erp).(condField).(grp){sj}.trans_counts = Transition_matrix;
                    Out.(erp).(condField).(grp){sj}.M_segments = M_segments;
                    Out.(erp).(condField).(grp){sj}.N_transitions = N_transitions;
                    Out.(erp).(condField).(grp){sj}.C_segments = C_segments;

                    for s = 1: 7
                        Out.(erp).(condField).(grp){sj}.(OccurrenceNames{s}) = C_segments(s) / max(1, M_segments);
                    end

                    Out.(erp).(condField).(grp){sj}.transition_count_matrix = Transition_matrix / max(1, N_transitions);
                    Out.(erp).(condField).(grp){sj}.transition_duration_matrix = Trans_dur;

                    clear T; T = Out.(erp).(condField).(grp){sj}.transition_duration_matrix; % For each subject, average over all trials the start and end duration for each transition

                    T_avg = cell(size(T));

                    for r = 1: size(T, 1)
                        for c = 1: size(T, 2)
                            clear xxx; xxx = T{r, c};

                            if isnumeric(xxx) && isscalar(xxx) && isnan(xxx)
                                T_avg{r, c} = [nan nan];

                            elseif isnumeric(xxx)
                                T_avg{r, c} = mean(xxx, 1, 'omitnan');

                            else
                                T_avg{r, c} = [nan nan];
                            end
                        end
                    end

                    Out.(erp).(condField).(grp){sj}.transition_duration_matrix = T_avg;

                end
            end
        end
    end

    [srcIdx, tgtIdx] = find(~eye(7)); % 42 directed pairs
    GLM_OPTS = statset('MaxIter', 5000, 'TolX', 1e-10, 'TolFun', 1e-10);

    LOGE_CLAMP = 12;
    MIN_OBS = 3;
    MIN_EXP = 3;
    Alpha = 0.05;

    %% ======================================================   %
    %    Q1. Deviations from independence within each (ERP×Group×Condition)  %
    %           Per-ERP BH-FDR post-pass pooling across all edges in the ERP.      %
    %    ======================================================= %

    TransStats = struct();

    % Pool Q1 p-values separately for each ERP

    poolQ1 = repmat(struct('p', [], 'map', []), numel(ERPs), 1);

    % map rows: [g c s t] for that ERP index (implicit by poolQ1(e))

    clear ALL_Mean; load([saveDir '\Metrics_ALL.mat']); % Load the dataset containing the three static measures in each condition in each subject in each group

    for e = 1: numel(ERPs)
        erp = ERPs{e};

        for g = 1: numel(Groups)
            grp = Groups{g};

            % z-score age within group (intercept interpretable at group-mean age)

            Ages_raw = AgeMap.(grp)(:);
            muA = mean(Ages_raw, 'omitnan');
            sdA = std(Ages_raw, 'omitnan');

            if sdA <= eps
                AgesZ = zeros(size(Ages_raw));
            else
                AgesZ = (Ages_raw - muA) ./ sdA;
            end

            for c = 1: numel(Conditions)
                condName = Conditions{c};
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

                RowsAll = table(); % y, e, logE, AgeZ, src, tgt

                if ~isempty(cellOut)
                    for sj = 1: numel(cellOut)

                        if sj > numel(AgesZ)
                            continue;
                        end

                        S = cellOut{sj};

                        if isempty(S) || ~isfield(S, 'trans_counts') || isempty(S.trans_counts)
                            continue;
                        end

                        if S.M_segments <= 0 || S.N_transitions <= 0
                            continue;
                        end

                        P = (S.C_segments(:) / S.M_segments);

                        Ssq = sum(P .^ 2);
                        den = 1 - Ssq;
                        kK = numel(S.C_segments);

                        if den <= eps
                            continue;
                        end

                        gamma = 1 / den;

                        E = gamma * S.N_transitions * (P * P.');
                        E(1: kK + 1: kK * kK) = NaN;

                        yvec = S.trans_counts(sub2ind([kK, kK], srcIdx, tgtIdx));
                        Dur_src_tgt = cell2mat(S.transition_duration_matrix(sub2ind([kK, kK], srcIdx, tgtIdx)));
                        Dur_src = Dur_src_tgt(:, 1);
                        Dur_trg = Dur_src_tgt(:, 2);
                        evec = E(sub2ind([kK, kK], srcIdx, tgtIdx));
                        keep = isfinite(evec) & evec > 0 & ~isnan(yvec);

                        if ~any(keep)
                            continue;
                        end

                        eKeep = evec(keep);
                        logE = log(eKeep);
                        logE = max(min(logE, LOGE_CLAMP), -LOGE_CLAMP);

                        clear R

                        R = table(yvec(keep), eKeep, logE, repmat(AgesZ(sj), nnz(keep), 1), srcIdx(keep), tgtIdx(keep), Dur_src(keep, :), Dur_trg(keep, :), 'VariableNames', {'y', 'e', 'logE', 'AgeZ', 'src', 'tgt', 'Dur_src', 'Dur_trg'});
                        RowsAll = [RowsAll; R];
                    end
                end

                for k = 1: numel(srcIdx)

                    s = srcIdx(k);
                    t = tgtIdx(k);
                    R = RowsAll(RowsAll.src == s & RowsAll.tgt == t, :);

                    if isempty(R)
                        continue;
                    end

                    if sum(R.y) < MIN_OBS || sum(R.e) < MIN_EXP
                        continue;
                    end

                    nRows(s, t) = height(R);
                    ObsOverExp(s, t) = sum(R.y) / sum(R.e);

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

                        if isfinite(pval(s, t))
                            poolQ1(e).p(end + 1, 1) = pval(s, t);
                            poolQ1(e).map(end + 1, :) = [g c s t];
                        end
                    catch
                    end
                end

                % Store raw results; per-ERP BH-FDR post-pass fills q, sig, dir, edges

                TransStats.(erp).(grp).(condName).IRR = IRR;
                TransStats.(erp).(grp).(condName).beta0 = beta0;
                TransStats.(erp).(grp).(condName).se = se;
                TransStats.(erp).(grp).(condName).p = pval;
                TransStats.(erp).(grp).(condName).IRR_CIlo = ciLo;
                TransStats.(erp).(grp).(condName).IRR_CIhi = ciHi;
                TransStats.(erp).(grp).(condName).nRows = nRows;
                TransStats.(erp).(grp).(condName).ObsOverExp = ObsOverExp;

                TransStats.(erp).(grp).(condName).pFDR = NaN(7);
                TransStats.(erp).(grp).(condName).sig = false(7);
                TransStats.(erp).(grp).(condName).dir = zeros(7);

                edges = table([], [], [], [], [], [], [], [], [], [], 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'ObsOverExp', 'p', 'pFDR', 'Dur_src', 'Dur_trg'});
                edges.srcLabel = strings(0, 1);
                edges.tgtLabel = strings(0, 1);
                edges.Direction = strings(0, 1);
                TransStats.(erp).(grp).(condName).edges = edges;

                TransStats.(erp).(grp).(condName).mFDR = NaN;
            end
        end
    end

    %% -------------------------- Q1: Per-ERP BH-FDR post-pass -------------------------- %%

    for e = 1: numel(ERPs)
        erp = ERPs{e};

        pPool = poolQ1(e).p;
        mFDR = numel(pPool);

        if mFDR <= 0
            continue;
        end

        pFDR_Pool = bh_adjust(pPool);

        % Assign edge-level pFDR back into each panel matrix

        for ii = 1: mFDR
            g = poolQ1(e).map(ii, 1);
            c = poolQ1(e).map(ii, 2);
            s = poolQ1(e).map(ii, 3);
            t = poolQ1(e).map(ii, 4);

            grp = Groups{g};
            condName = Conditions{c};

            TransStats.(erp).(grp).(condName).pFDR(s, t) = pFDR_Pool(ii);
        end

        % Build sig masks and edge lists per panel

        for g = 1: numel(Groups)
            grp = Groups{g};

            for c = 1: numel(Conditions)
                condName = Conditions{c};

                IRR = TransStats.(erp).(grp).(condName).IRR;
                ciLo = TransStats.(erp).(grp).(condName).IRR_CIlo;
                ciHi = TransStats.(erp).(grp).(condName).IRR_CIhi;
                pMat = TransStats.(erp).(grp).(condName).p;
                pFDRMat = TransStats.(erp).(grp).(condName).pFDR;
                ObsOverExp = TransStats.(erp).(grp).(condName).ObsOverExp;

                sigMask = isfinite(pFDRMat) & (pFDRMat <= 0.05);

                dirMatrix = zeros(7);
                dirMatrix(sigMask & IRR > 1) = 1;
                dirMatrix(sigMask & IRR < 1) = -1;

                [sSig, tSig] = find(sigMask);

                clear source_dur target_dur

                for pop = 1: length(sSig)
                    source_dur(pop) = nanmean(table2array(RowsAll(RowsAll.src == sSig(pop) & RowsAll.tgt == tSig(pop), 7)));
                    target_dur(pop) = nanmean(table2array(RowsAll(RowsAll.src == sSig(pop) & RowsAll.tgt == tSig(pop), 8)));
                end

                if ~isempty(sSig)
                    irr_vec = IRR(sub2ind([7, 7], sSig, tSig));
                    lo_vec = ciLo(sub2ind([7, 7], sSig, tSig));
                    hi_vec = ciHi(sub2ind([7, 7], sSig, tSig));
                    o2e_vec = ObsOverExp(sub2ind([7, 7], sSig, tSig));
                    p_vec = pMat(sub2ind([7, 7], sSig, tSig));
                    pFDR_vec = pFDRMat(sub2ind([7, 7], sSig, tSig));

                    edges = table(sSig, tSig, irr_vec, lo_vec, hi_vec, o2e_vec, p_vec, pFDR_vec, source_dur', target_dur', 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'ObsOverExp', 'p', 'pFDR', 'source_dur', 'target_dur'});
                    edges.srcLabel = reshape(string(MicrostateLabels(edges.src)), [], 1);
                    edges.tgtLabel = reshape(string(MicrostateLabels(edges.tgt)), [], 1);

                    dlab = strings(numel(sSig), 1);
                    dlab(irr_vec > 1) = "Above expected";
                    dlab(irr_vec < 1) = "Below expected";
                    edges.Direction = dlab;
                else
                    edges = table([], [], [], [], [], [], [], [], [], [], 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'ObsOverExp', 'p', 'pFDR', 'source_dur', 'target_dur'});
                    edges.srcLabel = strings(0, 1);
                    edges.tgtLabel = strings(0, 1);
                    edges.Direction = strings(0, 1);
                    edges.source_dur = strings(0, 1);
                    edges.target_dur = strings(0, 1);
                end

                TransStats.(erp).(grp).(condName).sig = sigMask;
                TransStats.(erp).(grp).(condName).dir = dirMatrix;
                TransStats.(erp).(grp).(condName).edges = edges;
                TransStats.(erp).(grp).(condName).mFDR = mFDR;
            end
        end
    end

    if ~Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
        save(fullfile(saveDir, ['Q1_base.mat']), 'TransStats');
    elseif ~Peaks_only_sensitivity_test && Smootness_penalty_sensitivity_test
        save(fullfile([saveDir 'Transition\Q1 and Q2'], ['Q1_sensitivity_smoothness.mat']), 'TransStats');
    elseif Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
        save(fullfile([saveDir 'Transition\Q1 and Q2'], ['Q1_sensitivity_peaks.mat']), 'TransStats');
    end

    %% ================================================ %%
    % Q2. Anchored contrasts (Target vs HC at the same age)
    %        Per-ERP BH-FDR is done inside buildAnchoredContrasts.m
    %% ================================================ %%

    CompStats = buildAnchoredContrasts(Out, AgeMap, ERPs, Conditions, Groups, 'RefGroup', 'HC', 'TargetGroups', {}, 'Alpha', 0.05, 'MaxIter', 5000, 'MinObs', MIN_OBS, 'MinExp', MIN_EXP, 'LogEClamp', LOGE_CLAMP);

    if ~Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
        save(fullfile(saveDir, ['Q2_base.mat']), 'CompStats');
    elseif ~Peaks_only_sensitivity_test && Smootness_penalty_sensitivity_test
        save(fullfile([saveDir 'Transition\Q1 and Q2'], ['Q2_sensitivity_smoothness.mat']), 'CompStats');
    elseif Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
        save(fullfile([saveDir 'Transition\Q1 and Q2'], ['Q2_sensitivity_peaks.mat']), 'CompStats');
    end

    if ~Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
        write_results_xlsx(saveDir, ERPs, Conditions, Groups, TransStats, CompStats);
    else
        write_results_xlsx([saveDir 'Transition\Q1 and Q2'], ERPs, Conditions, Groups, TransStats, CompStats);
    end

    % Compare baseline pipeline with smoothness panalty pipeline

    Conds = {'Negative', 'Neutral', 'Positive'};

    clear CompStatsTransStats
    if ~Peaks_only_sensitivity_test && Smootness_penalty_sensitivity_test
        load('D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\Q1_base.mat', 'TransStats'); Base = TransStats; clear TransStats
        load(fullfile(saveDir, ['\Transition\Q1 and Q2\Q1_sensitivity_smoothness.mat']), 'TransStats'); Sens = TransStats; clear TransStats
        [overall_Q1, byERP_Q1] = summarize_q1_robustness(Base, Sens, ERPs, Groups, Conds);

        load('D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\Q2_base.mat', 'CompStats'); Base = CompStats; clear CompStats
        load(fullfile(saveDir, ['\Transition\Q1 and Q2\Q2_sensitivity_smoothness.mat']), 'CompStats'); Sens = CompStats; clear CompStats
        [overall_Q2, byERP_Q2] = summarize_q2_robustness(Base, Sens, ERPs, Conds, Groups);

    elseif Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
        load('D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\Q1_base.mat', 'TransStats'); Base = TransStats; clear TransStats
        load(fullfile(saveDir, ['\Transition\Q1 and Q2\Q1_sensitivity_peaks.mat']), 'TransStats'); Sens = TransStats; clear TransStats
        [overall_Q1, byERP_Q1] = summarize_q1_robustness(Base, Sens, ERPs, Groups, Conds);

        load('D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\Q2_base.mat', 'CompStats'); Base = CompStats; clear CompStats
        load(fullfile(saveDir, ['\Transition\Q1 and Q2\Q2_sensitivity_peaks.mat']), 'CompStats'); Sens = CompStats; clear CompStats
        [overall_Q2, byERP_Q2] = summarize_q2_robustness(Base, Sens, ERPs, Conds, Groups);
    end

    %% ----------------------------- PLOTTING ------------------------------- %%

    plotting = 1;

    if plotting

        [WidthDomain, DurDomain] = computeGlobalWidthDomain(TransStats, CompStats, ERPs, Conditions, Groups);

        % Node radius and label

        NodeRadius = 0.2;
        NodeLabels = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};

        % Arrow width, color, starting and ending points

        MinArrowWidth = 0.3;
        MaxArrowWidth = 12.0;

        GreyColor = [0.90 0.90 0.90];
        GreyLineWidth = 1.0;
        AboveColor = [0.85 0.10 0.10];
        BelowColor = [0.10 0.35 0.95];

        % Bubbles radii, color, and edge

        DurRadiusMin = 0.03;
        DurRadiusMax = 0.13;
        DurFaceColor = [0.78 0.78 0.78];
        DurEdgeColor = [0 0 0];
        DurEdgeWidth = 0.5;

        Arguments = {'NodeLabels', NodeLabels, 'NodeRadius', NodeRadius, 'WidthDomain', WidthDomain, 'MinArrowWidth', MinArrowWidth, 'MaxArrowWidth', MaxArrowWidth, 'NodeLabels', NodeLabels, 'GreyColor', GreyColor, 'GreyLineWidth', GreyLineWidth, 'AboveColor', AboveColor, 'BelowColor', BelowColor, 'DurFaceColor', DurFaceColor, 'DurEdgeColor', DurEdgeColor, 'DurEdgeWidth', DurEdgeWidth, 'DurRadiusMin', DurRadiusMin, 'DurRadiusMax', DurRadiusMax, 'DurDomain', DurDomain};

        % ===============================================================================

        for e = 1: numel(ERPs)

            erp = ERPs{e};

            for g = 1: numel(Groups)
                grp = Groups{g};
                drawERPGridFigures(TransStats, erp, grp, Conditions, Arguments{:});
            end

            GroupContrastSlides(CompStats, erp, Conditions, Groups, Arguments{:});
        end
    end

end

if Link_trial_RT

    %% ===========================================================
    %  Link_trial_RT  — Reduced confirmatory family (UNCORRECTED)
    %  Trial-level log RT ~ syntax with within-subject×condition demeaning (Mundlak).
    %
    %  Primary tests (reported uncorrected):
    %   • N200 Anchor-2: stronger slope in valenced (Neg / Pos) vs Neutral?
    %   • P300 Anchor-2: stronger slope in valenced (Neg / Pos) vs Neutral?
    %
    %  Negative control (reported uncorrected): LPP backbone.
    %  Saves: RT_trial_confirmatory_*.csv
    %  ============================================================

    %% -----------------------------------
    % Load demographics and RT
    %% -----------------------------------

    load(fullfile(inputDir, 'France_data.mat'), 'France');

    %% ---------------------
    % Load outputStats
    %% ---------------------

    S = load(fullfile(saveDir, 'outputStats.mat'));
    outputStats = S.outputStats;

    anchors = findAnchorsBackbone(outputStats, ERP.window); % HC topology → anchors and LPP backbone (counts-based)
    subj_id_map = buildSubjectIdMap(France); % Subject-ID map (order alignment)

    %% -------------------------
    % 1) Per-trial RT table
    %% -------------------------

    rt_trials = buildRTTrialsTable(France); % subject, condition, trial, rt_ms, log_rt, age

    %% -----------------------------------------------------
    % 2) Per-trial motif predictors from MSClass
    %% -----------------------------------------------------

    trial_feats = computeTrialPredictors(outputStats, anchors, ERP.window, subj_id_map);

    %% -------------------------------
    % 3) Join and prepare data
    %% -------------------------------

    trial_data = innerjoin(rt_trials, trial_feats, 'Keys', {'subject', 'condition', 'trial'});
    trial_data.condition = categorical(trial_data.condition, {'Neutral', 'Negative', 'Positive'});

    % Ensure subject is categorical for fitlme grouping

    trial_data.subject = categorical(trial_data.subject);

    %% --------------------------------------------------
    % Standardize predictors (z) over all trials
    %% --------------------------------------------------

    pred_cols = {'N200_a1_rate', 'N200_a2_rate', 'P300_a1_rate', 'P300_a2_rate', 'LPP_backbone_rate'};

    for k = 1: numel(pred_cols)

        zname = ['z_' pred_cols{k}];
        mu = mean(trial_data.(pred_cols{k}), 'omitnan');
        sd = std(trial_data.(pred_cols{k}), 'omitnan');

        if ~isfinite(sd) || sd <= eps
            trial_data.(zname) = zeros(height(trial_data), 1);
        else
            trial_data.(zname) = (trial_data.(pred_cols{k}) - mu) ./ sd;
        end
    end

    %% ----------------------------------------------------
    % Standardize age across all included trials
    %% ----------------------------------------------------

    muA = mean(double(trial_data.age), 'omitnan');
    sdA = std(double(trial_data.age), 'omitnan');

    if ~isfinite(sdA) || sdA <= eps
        trial_data.age_z = zeros(height(trial_data), 1);
    else
        trial_data.age_z = (double(trial_data.age) - muA) ./ sdA;
    end

    %% ---------------------------------------------------------------
    % Standardize trial index within Subject × Condition
    %% ---------------------------------------------------------------

    [G_tc, ~] = findgroups(trial_data.subject, trial_data.condition);
    muT = splitapply(@mean, trial_data.trial, G_tc);
    sdT = splitapply(@std, trial_data.trial, G_tc);

    trial_data.trial_z = zeros(height(trial_data), 1);

    for r = 1: height(trial_data)

        g = G_tc(r);

        if ~isfinite(sdT(g)) || sdT(g) <= eps
            trial_data.trial_z(r) = 0;
        else
            trial_data.trial_z(r) = (trial_data.trial(r) - muT(g)) ./ sdT(g);
        end
    end

    %% ------------------------------------------------------
    % 4) Within subj×cond demeaning (Mundlak)
    %% ------------------------------------------------------

    pred_z = {'z_N200_a1_rate', 'z_N200_a2_rate', 'z_P300_a1_rate', 'z_P300_a2_rate', 'z_LPP_backbone_rate'};
    trial_data = demeanWithinSubjCond(trial_data, pred_z); % adds w_* and m_*

    %% -------------------------------------------------------------
    % 5) Valence contrast (Neutral = 0; Neg / Pos = 1)
    %% -------------------------------------------------------------

    trial_data.c_val = double(trial_data.condition ~= 'Neutral');

    %% ------------------------------------------------------
    % 6) Keep complete rows (required columns)
    %% ------------------------------------------------------

    need_base = {'z_N200_a2_rate', 'z_P300_a2_rate', 'z_LPP_backbone_rate'};
    use = isfinite(trial_data.log_rt) & isfinite(trial_data.age_z) & isfinite(trial_data.trial_z) & isfinite(trial_data.c_val);

    for k = 1: numel(need_base)
        use = use & isfinite(trial_data.(['w_' need_base{k}])) & isfinite(trial_data.(['m_' need_base{k}]));
    end

    trialW = trial_data(use, :);

    %% --------------------
    % Load outputStats
    %% --------------------

    S = load(fullfile(saveDir, 'outputStats.mat'));
    outputStats = S.outputStats;

    anchors = findAnchorsBackbone(outputStats, ERP.window); % HC topology → anchors and LPP backbone (counts-based)
    subj_id_map = buildSubjectIdMap(France); % Subject-ID map (order alignment)

    %% -----------------------------------------------------------------------------
    % 7) Fit the two Anchor-2 models (UNCORRECTED inference)
    %    (random slope for within-effect; safe fallback)
    %% -----------------------------------------------------------------------------

    R_N200 = fit_one_anchor2(trialW, 'z_N200_a2_rate', 'N200 Anchor-2 (valenced vs neutral)', anchors);
    R_P300 = fit_one_anchor2(trialW, 'z_P300_a2_rate', 'P300 Anchor-2 (valenced vs neutral)', anchors);

    confirm2 = [R_N200; R_P300];

    %% -------------------------------------------------
    % 8) Negative control (UNCORRECTED)
    %% -------------------------------------------------

    Rctl = fit_one_anchor2(trialW, 'z_LPP_backbone_rate', 'LPP backbone (control; uncorrected)', anchors);

    %% ---------------------------
    % 9) Save + print results
    %% ---------------------------

    outT = [confirm2; Rctl];

    if ~Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
        csvPath = fullfile(saveDir, 'RT_trial_confirmatory_base.csv');
        writetable(outT, csvPath);
        makeConfirmatory2Figure(saveDir);

    elseif Peaks_only_sensitivity_test && ~Smootness_penalty_sensitivity_test
        csvPath = fullfile([saveDir 'Transition\Q3'], 'RT_trial_confirmatory_only_peaks.csv');
        writetable(outT, csvPath);
        makeConfirmatory2Figure([saveDir 'Transition\Q3']);
        OUT_q3_smooth = summarize_q3_robustness('D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\RT_trial_confirmatory_base.csv', fullfile([saveDir 'Transition\Q3'], 'RT_trial_confirmatory_only_peaks.csv'), 'BaseLabel', "baseline", 'SensLabel', "only_peaks");

    elseif ~Peaks_only_sensitivity_test && Smootness_penalty_sensitivity_test
        csvPath = fullfile([saveDir 'Transition\Q3'], 'RT_trial_confirmatory_smoothness_penalty.csv');
        writetable(outT, csvPath);
        makeConfirmatory2Figure([saveDir 'Transition\Q3']);
        OUT_q3_smooth = summarize_q3_robustness('D:\Papers\2025\In preparation\XXX (task microstates)\Analysis\RT_trial_confirmatory_base.csv', fullfile([saveDir 'Transition\Q3'], 'RT_trial_confirmatory_smoothness_penalty.csv'), 'BaseLabel', "baseline", 'SensLabel', "smoothness_penalty");

    end

    disp('Reduced confirmatory family (UNCORRECTED): interaction w*c_val = slope increase in valenced vs neutral.');
    disp(outT(:, {'Predictor', 'Beta', 'CI_low', 'CI_high', 'pValue', 'Beta_NeutralSlope', 'pValue_NeutralSlope', 'Beta_ValencedSlope', 'pValue_ValencedSlope'}));

end

%% ====================================================
%% Fit one predictor: within × valence + Mundlak mean + covariates
%% ====================================================

function R = fit_one_anchor2(T, pred_base, pretty_name, anchors)

% Fit a single Anchor-2 (or control) model on within-demeaned data with valence contrast.
% Returns:
%  - Interaction (within × valence; confirmatory)
%  - Neutral slope (c_val = 0)
%  - Valenced slope (c_val = 1) computed as Neutral + Interaction
% Also annotates output with the HC-defined motif IDs used to build predictors.

if nargin < 4
    anchors = struct();
end

wcol = ['w_' pred_base]; % within-person deviation
mcol = ['m_' pred_base]; % Mundlak cluster mean (between-person)

% Keep complete rows for this predictor

use = isfinite(T.log_rt) & isfinite(T.age_z) & isfinite(T.trial_z) & isfinite(T.(wcol)) & isfinite(T.(mcol)) & isfinite(T.c_val);
TT = T(use, :);

% Fixed effects and random effects (try slope, else fallback)

fix = ['log_rt ~ ' wcol ' + c_val + ' wcol ':c_val + ' mcol ' + age_z + trial_z'];
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

% Compute coefficient CIs

CI = coefCI(lme, 'Alpha', 0.05);

% Confirmatory coefficient = interaction (within × valence)

term_inter = string([wcol ':c_val']);
idxI = strcmp(nam, term_inter);

% Neutral slope (c_val = 0) = main effect of within

idxN = strcmp(nam, string(wcol));

if ~any(idxI) || ~any(idxN)
    error('Could not find required fixed-effect terms in model coefficients: %s and %s.', wcol, term_inter);
end

% Valenced slope = Neutral + Interaction

b_neu = C.Estimate(idxN);
b_int = C.Estimate(idxI);
b_val = b_neu + b_int;

% SE for linear combination using covariance of fixed effects

CovB = lme.CoefficientCovariance;
ncoef = size(C, 1);

if size(CovB, 1) >= ncoef && size(CovB, 2) >= ncoef
    se_val = sqrt(CovB(idxN, idxN) + CovB(idxI, idxI) + 2 * CovB(idxN, idxI));
else
    se_val = NaN;
end

% p-value for linear combination (Neutral + Interaction)

L = zeros(1, ncoef);
L(idxN) = 1;
L(idxI) = 1;
p_val = coefTest(lme, L, 0);

% DF for CI on linear combination (conservative)

vnames = getVarNamesLocal(C);
hasDF = any(strcmpi(vnames, 'DF'));

if hasDF
    dfN = C.DF(idxN);
    dfI = C.DF(idxI);

    if all(isfinite([dfN dfI])) && all([dfN dfI] > 0)
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

% ----------------------
% Build output row
% ----------------------

R = table;

% Identifier columns

R.Predictor = string(pretty_name);
R.Formula = string(formula_try);
R.Term = term_inter;

% HC-defined motif IDs (saved into every row for provenance)

[R.N200_Anchor2_MapID, R.P300_Anchor2_MapID, R.LPP_Backbone_Map1, R.LPP_Backbone_Map2, R.LPP_Backbone_Pair] = anchors_to_columns(anchors);

% Interaction (confirmatory)

R.Beta = b_int;
R.CI_low = CI(idxI, 1);
R.CI_high = CI(idxI, 2);
R.pValue = C.pValue(idxI);

% Neutral slope

R.Beta_NeutralSlope = b_neu;
R.CI_low_NeutralSlope = CI(idxN, 1);
R.CI_high_NeutralSlope = CI(idxN, 2);
R.pValue_NeutralSlope = C.pValue(idxN);

% Valenced slope

R.Beta_ValencedSlope = b_val;
R.CI_low_ValencedSlope = ci_val(1);
R.CI_high_ValencedSlope = ci_val(2);
R.pValue_ValencedSlope = p_val;

end

function [n200_a2, p300_a2, lpp_i, lpp_j, lpp_pair] = anchors_to_columns(anchors)

% Default missing values

n200_a2 = NaN;
p300_a2 = NaN;
lpp_i = NaN;
lpp_j = NaN;
lpp_pair = "";

% N200 anchor-2

if isfield(anchors, 'N200') && numel(anchors.N200) >= 2
    n200_a2 = double(anchors.N200(2));
end

% P300 anchor-2

if isfield(anchors, 'P300') && numel(anchors.P300) >= 2
    p300_a2 = double(anchors.P300(2));
end

% LPP backbone pair

if isfield(anchors, 'LPP') && numel(anchors.LPP) >= 2
    lpp_i = double(anchors.LPP(1));
    lpp_j = double(anchors.LPP(2));
    lpp_pair = string(sprintf('%d <-> %d', lpp_i, lpp_j));
end

end

function vnames = getVarNamesLocal(C)

if istable(C)
    vnames = C.Properties.VariableNames;
elseif isa(C, 'dataset')
    vnames = get(C, 'VarNames');
else
    vnames = {};
end

end

%% ===========================================================
%% HC topology: two anchors per N200 / P300 and LPP reciprocal backbone
%% ===========================================================

function anchors = findAnchorsBackbone(outputStats, windows)

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

        X = subj_arr(s).MSClass;

        if isempty(X)
            continue;
        end

        counts_N200 = counts_N200 + transitionCounts(X, windows{1});
        counts_P300 = counts_P300 + transitionCounts(X, windows{2});
        counts_LPP = counts_LPP + transitionCounts(X, windows{3});
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

C = zeros(7, 7);
t1 = max(1, win(1));
t2 = min(size(MSClass, 1), win(2));

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

%% ==========================
%% Subject number alignment map
%% ==========================

function subj_id_map = buildSubjectIdMap(France)

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

%% ============
%% RT trials table
%% ============

function T = buildRTTrialsTable(France)

rows = {};

add_group_mood(France, 'BP_I', 'Depressed');
add_group_mood(France, 'BP_I', 'Euthymic');
add_group_mood(France, 'BP_II', 'Depressed');
add_group_mood(France, 'BP_II', 'Euthymic');
add_group_nomood(France, 'HC');
add_group_nomood(France, 'Siblings');

if isempty(rows)
    T = table('Size', [0 6], 'VariableTypes', {'double', 'string', 'double', 'double', 'double', 'double'}, 'VariableNames', {'subject', 'condition', 'trial', 'rt_ms', 'log_rt', 'age'});
else
    T = vertcat(rows{:});
end

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

        if ~isfield(s, cond) || ~isfield(s.(cond), 'Correct_resp_RT')
            return;
        end

        v = s.(cond).Correct_resp_RT;
        v = v(:).';

        for t = 1: numel(v)

            rt = v(t);

            if ~isfinite(rt)
                continue;
            end

            R = table;
            R.subject = double(s.number);
            R.condition = string(cond);
            R.trial = double(t);
            R.rt_ms = double(rt);
            R.log_rt = log(double(rt));
            R.age = double(s.age);

            rows{end + 1, 1} = R;
        end
    end

end

%% ================================
%% Per-trial motif predictors from MSClass
%% ================================

function T = computeTrialPredictors(outputStats, anchors, windows, subj_id_map)

rows = {};
Conditions = {'negative', 'neutral', 'positive'};
groups = fieldnames(outputStats);

for g = 1: numel(groups)

    grp = groups{g};

    if ~isfield(subj_id_map, grp)
        continue;
    end

    ids_vec = subj_id_map.(grp);

    for c = 1: numel(Conditions)

        cname = Conditions{c};

        if ~isfield(outputStats.(grp), cname)
            continue;
        end

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

if isempty(rows)
    T = table('Size', [0 8], 'VariableTypes', {'double', 'string', 'double', 'double', 'double', 'double', 'double', 'double'}, 'VariableNames', {'subject', 'condition', 'trial', 'N200_a1_rate', 'N200_a2_rate', 'P300_a1_rate', 'P300_a2_rate', 'LPP_backbone_rate'});
else
    T = vertcat(rows{:});
end

    function add_subject(rec, subjnum, cname)

        if ~isfinite(subjnum)
            return;
        end

        X = rec.MSClass;

        if isempty(X)
            return;
        end

        nT = size(X, 2);

        L.N200 = (windows{1}(2) - windows{1}(1) + 1) / 1000;
        L.P300 = (windows{2}(2) - windows{2}(1) + 1) / 1000;
        L.LPP = (windows{3}(2) - windows{3}(1) + 1) / 1000;

        L.N200 = max(1e-9, L.N200);
        L.P300 = max(1e-9, L.P300);
        L.LPP = max(1e-9, L.LPP);

        for t = 1: nT

            Cn = transitionCounts_single(X(windows{1}(1): windows{1}(2), t));
            a1 = anchors.N200(1);
            a2 = anchors.N200(2);
            n200_a1 = (sum(Cn(a1, :)) - Cn(a1, a1)) / L.N200;
            n200_a2 = (sum(Cn(a2, :)) - Cn(a2, a2)) / L.N200;

            Cp = transitionCounts_single(X(windows{2}(1): windows{2}(2), t));
            p1 = anchors.P300(1);
            p2 = anchors.P300(2);
            p300_a1 = (sum(Cp(p1, :)) - Cp(p1, p1)) / L.P300;
            p300_a2 = (sum(Cp(p2, :)) - Cp(p2, p2)) / L.P300;

            Cl = transitionCounts_single(X(windows{3}(1): windows{3}(2), t));
            i = anchors.LPP(1);
            j = anchors.LPP(2);
            lpp_bb = 0.5 * ((Cl(i, j) / L.LPP) + (Cl(j, i) / L.LPP));

            R = table;
            R.subject = double(subjnum);
            R.condition = string(upperFirst(cname));
            R.trial = double(t);
            R.N200_a1_rate = double(n200_a1);
            R.N200_a2_rate = double(n200_a2);
            R.P300_a1_rate = double(p300_a1);
            R.P300_a2_rate = double(p300_a2);
            R.LPP_backbone_rate = double(lpp_bb);

            rows{end + 1, 1} = R;
        end

    end

end

function C = transitionCounts_single(lbl)

C = zeros(7, 7);
l = double(lbl(:));
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

function s = upperFirst(sin)

s = lower(string(sin));

if strlength(s) == 0
    return;
end

ch = char(s);
ch(1) = upper(ch(1));
s = string(ch);

end

%% =======================================
%% Mundlak demeaning within Subject × Condition
%% =======================================

function T = demeanWithinSubjCond(T, pred_cols)

[G, ~] = findgroups(T.subject, T.condition);

for k = 1: numel(pred_cols)

    col = pred_cols{k};
    mu_gc = splitapply(@nanmean, T.(col), G);
    T.(['m_' col]) = mu_gc(G);
    T.(['w_' col]) = T.(col) - T.(['m_' col]);
end

end

%% ===================================
%% Confirmatory forest plot (UNCORRECTED)
%% ===================================

function outPath = makeConfirmatory2Figure(outputDir)

csv1 = fullfile(outputDir, 'RT_trial_confirmatory_base.csv');
csv2 = fullfile(outputDir, 'RT_trial_confirmatory_only_peaks.csv');
csv3 = fullfile(outputDir, 'RT_trial_confirmatory_smoothness_penalty.csv');

if exist(csv1, 'file')
    csvPath = csv1;
elseif exist(csv2, 'file')
    csvPath = csv2;
elseif exist(csv3, 'file')
    csvPath = csv3;
else
    error('Could not find RT_trial_confirmatory_*.csv in %s', outputDir);
end

T = readtable(csvPath);

keep = contains(string(T.Predictor), 'Anchor-2', 'IgnoreCase', true);
TT = T(keep, :);

if isempty(TT)
    error('No Anchor-2 rows found in %s', csvPath);
end

rowOrder = zeros(height(TT), 1);

for i = 1: height(TT)
    if contains(string(TT.Predictor(i)), 'P300', 'IgnoreCase', true)
        rowOrder(i) = 1;
    else
        rowOrder(i) = 2;
    end
end

[~, ix] = sort(rowOrder);
TT = TT(ix, :);

n = height(TT);
yBase = (1: n).';

bN = TT.Beta_NeutralSlope;
ciN = [TT.CI_low_NeutralSlope TT.CI_high_NeutralSlope];
bV = TT.Beta_ValencedSlope;
ciV = [TT.CI_low_ValencedSlope TT.CI_high_ValencedSlope];

colNeutral = [0.60 0.60 0.60];
colVal = [0.85 0.20 0.20];
mkN = 'o';
mkV = 's';
off = [-0.17, +0.17];

fig = figure('Color', 'w', 'Position', [120 120 780 360]);
hold on;

for i = 1: n

    yN = yBase(i) + off(1);
    yV = yBase(i) + off(2);

    plot(ciN(i, :), [yN yN], '-', 'Color', colNeutral, 'LineWidth', 2);
    plot(bN(i), yN, mkN, 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colNeutral);

    plot(ciV(i, :), [yV yV], '-', 'Color', colVal, 'LineWidth', 2);

    fillVal = 'w';

    if ismember('pValue', TT.Properties.VariableNames) && isfinite(TT.pValue(i)) && TT.pValue(i) < 0.05
        fillVal = colVal;
    end

    plot(bV(i), yV, mkV, 'MarkerSize', 7, 'MarkerFaceColor', fillVal, 'MarkerEdgeColor', colVal);
end

xline(0, ':', 'Color', [0.4 0.4 0.4]);

labs = strings(n, 1);

for i = 1: n
    if contains(string(TT.Predictor(i)), 'P300', 'IgnoreCase', true)
        labs(i) = "P300 Anchor";
    else
        labs(i) = "N200 Anchor";
    end
end

set(gca, 'YTick', yBase, 'YTickLabel', labs, 'YDir', 'normal', 'Box', 'off');
xlabel('Standardized \beta on log RT (trial-level, within-subject×condition centered)');
legend({'Neutral', 'Valenced (Neg / Pos)'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
title('Behavior–syntax effects (Neutral vs Valenced slopes)\newline(Filled = uncorrected p < 0.05 on interaction)');
grid on;

outPath = fullfile(outputDir, 'Fig_RT_Confirmatory.tif');
exportgraphics(fig, outPath, 'Resolution', 300, 'BackgroundColor', 'white');
close(fig);

fprintf('Saved %s\n', outPath);

end
