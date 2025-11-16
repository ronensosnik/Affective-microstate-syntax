clear all
close all hidden
warning off
clc

Groups = {'BP_I_Depressed', 'BP_I_Euthymic', 'BP_II_Depressed', 'BP_II_Euthymic', 'HC', 'Siblings'};

% What actions to perform

Identify_individ_microstate_maps = 0;
Identify_mean_microstate_maps_for_each_group = 0;
Identify_grand_mean_microstate_maps = 0;
Sort_grand_mean_maps = 0;
Plot_grand_mean_maps = 0;
Compare_classes_across_groups = 0;
Compare_grand_mean_Custo = 0;
Backfit_and_quantify_temp_dynam = 0;

Prepare_demogarphics_table = 0;

Analyze = 1;

%% Setting all parameters

test_backfitting_penalty = 0;

% Set path to directory containing group folders

inputDir = 'D:\Papers\2025\Submitted\Journal of Bipolar Disorders (resting state microstates)\Data';

if ~test_backfitting_penalty
    saveDir = 'D:\Papers\2025\Submitted\Journal of Bipolar Disorders (resting state microstates)\Analysis';
else
    saveDir = 'D:\';
end

groupFolders = dir(inputDir);
groupFolders = groupFolders([groupFolders.isdir] & ~ismember({groupFolders.name},{'.', '..'}));
nGroups = length(groupFolders);
groupNames = cell(1, nGroups);
condNames = cell(1, nGroups);
dataDirs = cell(1, nGroups);

% Set clustering parameters

ClustPar.UseAAHC = false;            % true = AAHC, false = kmeans
ClustPar.MinClasses = 3;                % minimum number of clusters to identify
ClustPar.MaxClasses = 15;             % maximum number of clusters to identify
ClustPar.Restarts = 100;                 % number of times kmeans algorithm is restarted (ignored if using AAHC)
ClustPar.Allow_early_stop = true;   % terminating when the best GEV in the most recent batch failed to exceed the running best by ≥0.01%
ClustPar.MaxMaps = inf;                 % maximum number of data samples to use to identify clusters
ClustPar.GFPPeaks = false;            % whether clustering should be limited to global field power peaks
ClustPar.IgnorePolarity = true;         % whether maps of inverted polarities should be considered part of the same cluster
ClustPar.Normalize = true;               % Set to false if using AAHC

% Set backfitting parameters

FitPar.Classes = 7;    % cluster solutions to use for backfitting
FitPar.PeakFit = 0;     % whether to backfit only on global field power peaks
FitPar.lambda = 0.3;  % smoothness penalty - ignored if FitPar.PeakFit = 1
FitPar.b = 30;             % smoothing window (ms) - ignored if FitPar.PeakFit = 1

if test_backfitting_penalty
    FitPar.lambda = 0.8;
end

for i=1: nGroups
    groupDir = fullfile(inputDir, groupFolders(i).name);
    groupNames{i} = groupFolders(i).name;
    condFolders = dir(groupDir);
    condFolders = condFolders(~matches({condFolders.name}, {'.', '..'}));
    Subj_num{i} = {condFolders.name};
    dataDirs{i} = cellfun(@(x) fullfile(groupDir, x), Subj_num{i}, 'UniformOutput', false);
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
    MSTemplate = eeg_store(MSTemplate,pop_loadset('filename',Templates(t).name,'filepath',templatepath));
end

global MSTEMPLATE;
MSTEMPLATE = MSTemplate;

GroupIdx = cell(1, nGroups);
lastGroupIdx = 1;

%% Identify_individ_microstate_maps

if Identify_individ_microstate_maps

    % Load datasets and update subject, group, and condition info

    for i = 1: nGroups

        for j = 1: numel(dataDirs{i})
            setFiles = dir(fullfile(dataDirs{i}{j}, '*.set'));
            setFilenames = {setFiles.name};

            % Load datasets

            fprintf('Loading datasets in group %s...\n', groupNames{i});
            EEG = pop_loadset('filename', setFilenames, 'filepath', dataDirs{i}{j});
            [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET);
            currGroupIdx = lastGroupIdx: numel(ALLEEG);

            % Update group and condition info for all sets

            fprintf('Updating group and condition information for group %s, condition %s...\n', groupNames{i}, Subj_num{i}{j});

            for k = 1: numel(currGroupIdx)
                [EEG, ALLEEG, CURRENTSET] = eeg_retrieve(ALLEEG, currGroupIdx(k));
                filename = EEG.filename(1:strfind(EEG.filename, '.') - 1);
                idx = strfind(filename, '_');

                if isempty(idx)
                    EEG.subject = filename;
                else
                    EEG.subject = filename(1: idx(2) - 1);
                end

                EEG.group = groupNames{i};
                [ALLEEG, EEG,CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
            end

            GroupIdx{i} = [GroupIdx{i} currGroupIdx];
            lastGroupIdx = numel(ALLEEG) + 1;
        end
    end

    AllSubjects = 1: numel(ALLEEG);

    disp('Identifying microstates for all sets...');

    if (ClustPar.MinClasses == 7 && ClustPar.MaxClasses == 7)
        [EEG, CURRENTSET] = pop_FindMSMaps(ALLEEG, AllSubjects, 'ClustPar', ClustPar);
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        save ([saveDir '\Individual_micro_states.mat'], 'ALLEEG', 'EEG', 'CURRENTSET', '-v7.3');

    elseif (ClustPar.MinClasses == 3 && ClustPar.MaxClasses == 15)
        [ALLEEG, GEV, MeanGEV, StdGEV, Ks, groupOrder] = findMSMapsAndGEV_bySubject_parallel(ALLEEG, ClustPar, MSTEMPLATE, string(pluginpath));
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
        save ([saveDir '\Individual_micro_states_3_15.mat'], 'ALLEEG', 'EEG', 'GEV', 'MeanGEV', 'StdGEV', 'Ks', 'groupOrder', 'CURRENTSET', '-v7.3');
        Compute_and_plot_individ_total_GEV(GEV)
    end

end

%% Identify_mean_microstate_maps_for_each_group

if Identify_mean_microstate_maps_for_each_group

    load('D:\Papers\2025\Submitted\Journal of Bipolar Disorders (resting state microstates)\Analysis\Individual_micro_states.mat');

    GroupMeanIdx = [];

    for i = 1: nGroups

        fprintf('Identifying mean maps for group %s\n', groupNames{i});

        EEG = pop_CombMSMaps(ALLEEG, GroupIdx{i}, 'MeanName', ['Mean_' groupNames{i}], 'IgnorePolarity', ClustPar.IgnorePolarity);
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');
        GroupMeanIdx = [GroupMeanIdx CURRENTSET];
    end

    save ([saveDir '\Individual_micro_states.mat'], 'ALLEEG', 'EEG', 'CURRENTSET', 'GroupMeanIdx', '-v7.3');

end

%% Identify_grand_mean_microstate_maps

if Identify_grand_mean_microstate_maps

    fprintf('Identifying grand mean maps');

    EEG = pop_CombMSMaps(ALLEEG, GroupMeanIdx, 'MeanName', 'Grand_mean', 'IgnorePolarity', ClustPar.IgnorePolarity);
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'gui', 'off');

    save ([saveDir '\Individual_micro_states.mat'], 'ALLEEG', 'EEG', 'CURRENTSET', '-v7.3');

end

%% Sort grand mean maps by the 2017 Custo maps and then sort all individual maps by the grand mean map

if Sort_grand_mean_maps

    % Sort the mean maps by the specified published template(s)

    disp(['Grand mean (across groups) mean (across seven classes) shared variance: ' num2str(round(mean(ALLEEG(CURRENTSET).msinfo.MSMaps(7).SharedVar) * 10000) / 100)]);

    % First, sort the 7-class solution of the grand mean by the 2017 Custo maps.

    [ALLEEG, EEG, CURRENTSET] = pop_SortMSMaps(ALLEEG, CURRENTSET, 'TemplateSet', 'Custo2017', 'Classes', 7, 'IgnorePolarity', 1);

    % Next, sort the seven classes of each individual template in each group by the seven classes of the grand mean

    [ALLEEG, EEG, CURRENTSET] = pop_SortMSMaps(ALLEEG, 1: CURRENTSET - 1, 'TemplateSet', CURRENTSET, 'Classes', 7, 'IgnorePolarity', 1);

    save ([saveDir '\Individual_micro_states.mat'], 'ALLEEG', 'EEG', 'CURRENTSET', '-v7.3');

end

if Plot_grand_mean_maps

    % Plot the sorted, seven classes of the grand mean

    eeglab;
    load([saveDir '\Individual_micro_states.mat']);
    pop_ShowIndMSMaps(ALLEEG, 67: 73, 'Classes', 7);

end

%% Compare classes across groups

if Compare_classes_across_groups

    load([saveDir '\Individual_micro_states.mat']);

    % Compare group means

    sharedVarTable = pop_CompareMSMaps(ALLEEG, 'MeanSets', 67: 72, 'Classes', 7, 'gui', 0);
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

    disp(['A_microstate = ' num2str(A_microstate_mean) ' +- ' num2str(A_microstate_std)]);
    disp(['B_microstate = ' num2str(B_microstate_mean) ' +- ' num2str(B_microstate_std)]);
    disp(['C_microstate = ' num2str(C_microstate_mean) ' +- ' num2str(C_microstate_std)]);
    disp(['D_microstate = ' num2str(D_microstate_mean) ' +- ' num2str(D_microstate_std)]);
    disp(['E_microstate = ' num2str(E_microstate_mean) ' +- ' num2str(E_microstate_std)]);
    disp(['F_microstate = ' num2str(F_microstate_mean) ' +- ' num2str(F_microstate_std)]);
    disp(['G_microstate = ' num2str(G_microstate_mean) ' +- ' num2str(G_microstate_std)]);

    disp(' ');
    disp(' ');
    disp(' ');

    % Compare grand mean and Custo et al template

    sharedVarTable_GrandMean = pop_CompareMSMaps(ALLEEG, 'MeanSets', 73, 'PublishedSets', 'Custo2017', 'Classes', 7, 'Filename', [saveDir '\sharedvars.mat'], 'gui', 0);
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

    sharedVarTablee = sharedVarTable_GrandMean;

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

    figure(2); set(gcf, 'Color', [1 1 1]); heatmap(sharedVarTable_plot, 'CellLabelColor', 'none', 'MissingDataColor', [1 1 1], 'GridVisible', 'off', 'ColorBarVisible' ,'off', 'Title', 'Comparison. Grand mean vs. Custo template'); colormap('jet');
    colorbar;

end

if Backfit_and_quantify_temp_dynam

    %% Backfit using Custo2017 and quantify temporal dynamics

    Groups = {'BP_II_Depressed', 'BP_II_Euthymic', 'BP_I_Depressed', 'BP_I_Euthymic', 'HC', 'Siblings'};
    load([saveDir '\Individual_micro_states.mat']);

    disp('Backfitting and extracting temporal dynamics...');

    [EEG, CURRENTSET] = pop_FitMSMaps(ALLEEG, 1: 66, 'TemplateSet', 'Custo2017', 'FitPar', FitPar); % The function already kills microstates truncated by boundaries
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);

    index = 0;

    for pop = 1: 12
        index = index + 1;
        outputStats.BP_II_Depressed(index).MSClass = ALLEEG(pop).msinfo.MSStats(7).MSClass;
        outputStats.BP_II_Depressed(index).GFP = ALLEEG(pop).msinfo.MSStats(7).GFP;
        outputStats.BP_II_Depressed(index).MeanDuration = ALLEEG(pop).msinfo.MSStats(7).MeanDuration;
        outputStats.BP_II_Depressed(index).MeanOccurrence = ALLEEG(pop).msinfo.MSStats(7).MeanOccurrence;
        outputStats.BP_II_Depressed(index).Coverage = ALLEEG(pop).msinfo.MSStats(7).Coverage;
    end

    index = 0;

    for pop = 13: 24
        index = index + 1;
        outputStats.BP_II_Euthymic(index).MSClass = ALLEEG(pop).msinfo.MSStats(7).MSClass;
        outputStats.BP_II_Euthymic(index).GFP = ALLEEG(pop).msinfo.MSStats(7).GFP;
        outputStats.BP_II_Euthymic(index).MeanDuration = ALLEEG(pop).msinfo.MSStats(7).MeanDuration;
        outputStats.BP_II_Euthymic(index).MeanOccurrence = ALLEEG(pop).msinfo.MSStats(7).MeanOccurrence;
        outputStats.BP_II_Euthymic(index).Coverage = ALLEEG(pop).msinfo.MSStats(7).Coverage;
    end

    index = 0;

    for pop = 25: 33
        index = index + 1;
        outputStats.BP_I_Depressed(index).MSClass = ALLEEG(pop).msinfo.MSStats(7).MSClass;
        outputStats.BP_I_Depressed(index).GFP = ALLEEG(pop).msinfo.MSStats(7).GFP;
        outputStats.BP_I_Depressed(index).MeanDuration = ALLEEG(pop).msinfo.MSStats(7).MeanDuration;
        outputStats.BP_I_Depressed(index).MeanOccurrence = ALLEEG(pop).msinfo.MSStats(7).MeanOccurrence;
        outputStats.BP_I_Depressed(index).Coverage = ALLEEG(pop).msinfo.MSStats(7).Coverage;
    end

    index = 0;

    for pop = 34: 42
        index = index + 1;
        outputStats.BP_I_Euthymic(index).MSClass = ALLEEG(pop).msinfo.MSStats(7).MSClass;
        outputStats.BP_I_Euthymic(index).GFP = ALLEEG(pop).msinfo.MSStats(7).GFP;
        outputStats.BP_I_Euthymic(index).MeanDuration = ALLEEG(pop).msinfo.MSStats(7).MeanDuration;
        outputStats.BP_I_Euthymic(index).MeanOccurrence = ALLEEG(pop).msinfo.MSStats(7).MeanOccurrence;
        outputStats.BP_I_Euthymic(index).Coverage = ALLEEG(pop).msinfo.MSStats(7).Coverage;
    end

    index = 0;

    for pop = 43: 54
        index = index + 1;
        outputStats.HC(index).MSClass = ALLEEG(pop).msinfo.MSStats(7).MSClass;
        outputStats.HC(index).GFP = ALLEEG(pop).msinfo.MSStats(7).GFP;
        outputStats.HC(index).MeanDuration = ALLEEG(pop).msinfo.MSStats(7).MeanDuration;
        outputStats.HC(index).MeanOccurrence = ALLEEG(pop).msinfo.MSStats(7).MeanOccurrence;
        outputStats.HC(index).Coverage = ALLEEG(pop).msinfo.MSStats(7).Coverage;
    end

    index = 0;

    for pop = 55: 66
        index = index + 1;
        outputStats.Siblings(index).MSClass = ALLEEG(pop).msinfo.MSStats(7).MSClass;
        outputStats.Siblings(index).GFP = ALLEEG(pop).msinfo.MSStats(7).GFP;
        outputStats.Siblings(index).MeanDuration = ALLEEG(pop).msinfo.MSStats(7).MeanDuration;
        outputStats.Siblings(index).MeanOccurrence = ALLEEG(pop).msinfo.MSStats(7).MeanOccurrence;
        outputStats.Siblings(index).Coverage = ALLEEG(pop).msinfo.MSStats(7).Coverage;
    end

    save(fullfile(saveDir, ['outputStats.mat']), 'outputStats');

end

if Prepare_demogarphics_table

    for i = 1: nGroups

        for j = 1: numel(dataDirs{i})

            clear CRF_dir CRF_Files_name

            CRF_dir = dir(fullfile(dataDirs{i}{j}, '*.mat'));
            CRF_Files_name = CRF_dir.name;

            load([CRF_dir.folder '\' CRF_Files_name]);

            switch groupNames{i}

                case 'BP_I_Depressed'

                    France.BP_I.Depressed{j}.age = Subject.General_info.Age_at_examination;
                    France.BP_I.Depressed{j}.gender = Subject.General_info.sex;
                    France.BP_I.Depressed{j}.number = Subject.General_info.number;
                    France.BP_I.Depressed{j}.GAF = Subject.GAF.score;
                    France.BP_I.Depressed{j}.MADRS = sum([Subject.MADRS.Apparent_sadness Subject.MADRS.Reported_sadness Subject.MADRS.Inner_tension Subject.MADRS.Reduced_sleep Subject.MADRS.Reduced_appetite Subject.MADRS.Concentration_difficulties Subject.MADRS.Lassitude Subject.MADRS.Inability_to_feel Subject.MADRS.Pessimistic_thoughts Subject.MADRS.Suicidal_thoughts]);
                    France.BP_I.Depressed{j}.YMRS = sum([Subject.MARS.Elevated_mood Subject.MARS.Increased_motor_activity_energy Subject.MARS.Sexual_interest Subject.MARS.Sleep Subject.MARS.Irritability Subject.MARS.Speech Subject.MARS.Language_thought_disorder Subject.MARS.Content Subject.MARS.Disruptive_aggressive_behavior Subject.MARS.Appearance Subject.MARS.Insight]);

                case 'BP_I_Euthymic'

                    France.BP_I.Euthymic{j}.age = Subject.General_info.Age_at_examination;
                    France.BP_I.Euthymic{j}.gender = Subject.General_info.sex;
                    France.BP_I.Euthymic{j}.number = Subject.General_info.number;
                    France.BP_I.Euthymic{j}.GAF = Subject.GAF.score;
                    France.BP_I.Euthymic{j}.MADRS = sum([Subject.MADRS.Apparent_sadness Subject.MADRS.Reported_sadness Subject.MADRS.Inner_tension Subject.MADRS.Reduced_sleep Subject.MADRS.Reduced_appetite Subject.MADRS.Concentration_difficulties Subject.MADRS.Lassitude Subject.MADRS.Inability_to_feel Subject.MADRS.Pessimistic_thoughts Subject.MADRS.Suicidal_thoughts]);
                    France.BP_I.Euthymic{j}.YMRS = sum([Subject.MARS.Elevated_mood Subject.MARS.Increased_motor_activity_energy Subject.MARS.Sexual_interest Subject.MARS.Sleep Subject.MARS.Irritability Subject.MARS.Speech Subject.MARS.Language_thought_disorder Subject.MARS.Content Subject.MARS.Disruptive_aggressive_behavior Subject.MARS.Appearance Subject.MARS.Insight]);

                case 'BP_II_Depressed'

                    France.BP_II.Depressed{j}.age = Subject.General_info.Age_at_examination;
                    France.BP_II.Depressed{j}.gender = Subject.General_info.sex;
                    France.BP_II.Depressed{j}.number = Subject.General_info.number;
                    France.BP_II.Depressed{j}.GAF = Subject.GAF.score
                    France.BP_II.Depressed{j}.MADRS = sum([Subject.MADRS.Apparent_sadness Subject.MADRS.Reported_sadness Subject.MADRS.Inner_tension Subject.MADRS.Reduced_sleep Subject.MADRS.Reduced_appetite Subject.MADRS.Concentration_difficulties Subject.MADRS.Lassitude Subject.MADRS.Inability_to_feel Subject.MADRS.Pessimistic_thoughts Subject.MADRS.Suicidal_thoughts]);
                    France.BP_II.Depressed{j}.YMRS = sum([Subject.MARS.Elevated_mood Subject.MARS.Increased_motor_activity_energy Subject.MARS.Sexual_interest Subject.MARS.Sleep Subject.MARS.Irritability Subject.MARS.Speech Subject.MARS.Language_thought_disorder Subject.MARS.Content Subject.MARS.Disruptive_aggressive_behavior Subject.MARS.Appearance Subject.MARS.Insight]);

                case 'BP_II_Euthymic'

                    France.BP_II.Euthymic{j}.age = Subject.General_info.Age_at_examination;
                    France.BP_II.Euthymic{j}.gender = Subject.General_info.sex;
                    France.BP_II.Euthymic{j}.number = Subject.General_info.number;
                    France.BP_II.Euthymic{j}.GAF = Subject.GAF.score;
                    France.BP_II.Euthymic{j}.MADRS = sum([Subject.MADRS.Apparent_sadness Subject.MADRS.Reported_sadness Subject.MADRS.Inner_tension Subject.MADRS.Reduced_sleep Subject.MADRS.Reduced_appetite Subject.MADRS.Concentration_difficulties Subject.MADRS.Lassitude Subject.MADRS.Inability_to_feel Subject.MADRS.Pessimistic_thoughts Subject.MADRS.Suicidal_thoughts]);
                    France.BP_II.Euthymic{j}.YMRS = sum([Subject.MARS.Elevated_mood Subject.MARS.Increased_motor_activity_energy Subject.MARS.Sexual_interest Subject.MARS.Sleep Subject.MARS.Irritability Subject.MARS.Speech Subject.MARS.Language_thought_disorder Subject.MARS.Content Subject.MARS.Disruptive_aggressive_behavior Subject.MARS.Appearance Subject.MARS.Insight]);

                case 'HC'

                    France.HC{j}.age = Subject.General_info.Age_at_examination;
                    France.HC{j}.gender = Subject.General_info.sex;
                    France.HC{j}.number = Subject.General_info.number;
                    France.HC{j}.GAF = 0/0;
                    France.HC{j}.MADRS = 0/0;
                    France.HC{j}.YMRS = 0/0;

                case 'Siblings'

                    France.Siblings{j}.age = Subject.General_info.Age_at_examination;
                    France.Siblings{j}.gender = Subject.General_info.sex;
                    France.Siblings{j}.number = Subject.General_info.number;
                    France.Siblings{j}.GAF = 0/0;
                    France.Siblings{j}.MADRS = 0/0;
                    France.Siblings{j}.YMRS = 0/0;

            end

        end
    end

    save('D:\Papers\2025\Submitted\Journal of Bipolar Disorders (resting state microstates)\Data\Participants_demographics', 'France');
end

if Analyze

    Dataset_info = 0;
    Check_demog = 0;
    analyze_metrics = 0;
    test_transition_probability = 1;

    % Find how many channels and PCA were omitted

    if Dataset_info

        for i = 1: nGroups

            for j = 1: numel(dataDirs{i})

                setFiles = dir([dataDirs{i}{j} '\Not_interpolated\*.set']);
                setFilenames = {setFiles.name};

                clear EEG; EEG = pop_loadset('filename', setFilenames, 'filepath', [dataDirs{i}{j} '\Not_interpolated\']);

                try
                    num_elec.(groupNames{i}) = [num_elec.(groupNames{i}) EEG.nbchan];
                catch
                    num_elec.(groupNames{i}) = EEG.nbchan;
                end

                try
                    num_ICA.(groupNames{i}) = [num_ICA.(groupNames{i}) size(EEG.icaweights, 1)];
                catch
                    num_ICA.(groupNames{i}) = size(EEG.icaweights, 1);
                end

            end

            Removed_elec.(groupNames{i}) = 62 - num_elec.(groupNames{i});
            Removed_ICA_comp.(groupNames{i}) = num_elec.(groupNames{i}) - num_ICA.(groupNames{i});
            Removed_ICA_comp.(groupNames{i})(find(Removed_ICA_comp.(groupNames{i}) < 0)) = 0;
        end

        All_removed_elec = [Removed_elec.(groupNames{1}) Removed_elec.(groupNames{2}) Removed_elec.(groupNames{3}) Removed_elec.(groupNames{4}) Removed_elec.(groupNames{5}) Removed_elec.(groupNames{6})]
        All_removed_ICA_comp = [Removed_ICA_comp.(groupNames{1}) Removed_ICA_comp.(groupNames{2}) Removed_ICA_comp.(groupNames{3}) Removed_ICA_comp.(groupNames{4}) Removed_ICA_comp.(groupNames{5}) Removed_ICA_comp.(groupNames{6})]

        % Electrodes

        disp([num2str(100 - (sum(All_removed_elec == 0) / 66) * 100) ' of the datasets required any channel removal']);
        disp(['range across groups ' num2str(min(All_removed_elec)) ' - ' num2str(max(All_removed_elec))]);

        for pop = 1: 6
            disp(['median removed by group: '  num2str(median(Removed_elec.(groupNames{pop}))) ' in ' groupNames{pop}]);
        end

        disp(' ');
        disp(' ');

        % ICA components

        disp([num2str(100 - (sum(All_removed_ICA_comp == 0) / 66) * 100) ' of the datasets required any component removal']);
        disp(['range across groups ' num2str(min(All_removed_ICA_comp)) ' - ' num2str(max(All_removed_ICA_comp))]);

        for pop = 1: 6
            disp(['median removed by group: ' num2str(median(Removed_ICA_comp.(groupNames{pop}))) ' in ' groupNames{pop}]);
        end
    end

    if Check_demog

        % Check demographics (age, gender, GAF, MADRS, YMRS)

        load('D:\Papers\2025\Submitted\Journal of Bipolar Disorders (resting state microstates)\Data\Participants_demographics.mat');
        Check_demographics(France);

    end

    if analyze_metrics

        Sens_meds = 1; medSpec = 'poly';

        % We want to test for an effect of group on mean microstate duration, microstate occurrence and coverage
        % while controling for age and gender

        % First, we load the information of each subject (age, gender, and other scores)

        load('D:\Papers\2025\Submitted\Journal of Bipolar Disorders (resting state microstates)\Data\Participants_demographics.mat');

        % Then, we load the microstate metrics (stored in outputStats) of each subject in each of the six groups

        load([saveDir '\outputStats.mat']);

        Duration_A = []; Occurrence_A = []; Coverage_A = [];
        Duration_B = []; Occurrence_B = []; Coverage_B = [];
        Duration_C = []; Occurrence_C = []; Coverage_C = [];
        Duration_D = []; Occurrence_D = []; Coverage_D = [];
        Duration_E = []; Occurrence_E = []; Coverage_E = [];
        Duration_F = []; Occurrence_F = []; Coverage_F = [];
        Duration_G = []; Occurrence_G = []; Coverage_G = [];

        Group = {};
        Age = []; Gender = {}; SubjectID = [];
        GAF = []; MADRS = []; YMRS = [];
        AD = []; AP = []; MS = []; ANX = []; OTHER = []; poly_count = []; poly_ge2 = []; drug_count = []; meds_missing = [];

        % Now, pool the Duration, Occurrence and Coverage of each of the seven microstates for each subject in each group

        for subject = 1: size(outputStats.BP_I_Depressed, 2) % Start with BP_I_Depressed

            Duration_A = [Duration_A outputStats.BP_I_Depressed(subject).MeanDuration(1)];
            Occurrence_A = [Occurrence_A outputStats.BP_I_Depressed(subject).MeanOccurrence(1)];
            Coverage_A = [Coverage_A outputStats.BP_I_Depressed(subject).Coverage(1)];

            Duration_B = [Duration_B outputStats.BP_I_Depressed(subject).MeanDuration(2)];
            Occurrence_B = [Occurrence_B outputStats.BP_I_Depressed(subject).MeanOccurrence(2)];
            Coverage_B = [Coverage_B outputStats.BP_I_Depressed(subject).Coverage(2)];

            Duration_C = [Duration_C outputStats.BP_I_Depressed(subject).MeanDuration(3)];
            Occurrence_C = [Occurrence_C outputStats.BP_I_Depressed(subject).MeanOccurrence(3)];
            Coverage_C = [Coverage_C outputStats.BP_I_Depressed(subject).Coverage(3)];

            Duration_D = [Duration_D outputStats.BP_I_Depressed(subject).MeanDuration(4)];
            Occurrence_D = [Occurrence_D outputStats.BP_I_Depressed(subject).MeanOccurrence(4)];
            Coverage_D = [Coverage_D outputStats.BP_I_Depressed(subject).Coverage(4)];

            Duration_E = [Duration_E outputStats.BP_I_Depressed(subject).MeanDuration(5)];
            Occurrence_E = [Occurrence_E outputStats.BP_I_Depressed(subject).MeanOccurrence(5)];
            Coverage_E = [Coverage_E outputStats.BP_I_Depressed(subject).Coverage(5)];

            Duration_F = [Duration_F outputStats.BP_I_Depressed(subject).MeanDuration(6)];
            Occurrence_F = [Occurrence_F outputStats.BP_I_Depressed(subject).MeanOccurrence(6)];
            Coverage_F = [Coverage_F outputStats.BP_I_Depressed(subject).Coverage(6)];

            Duration_G = [Duration_G outputStats.BP_I_Depressed(subject).MeanDuration(7)];
            Occurrence_G = [Occurrence_G outputStats.BP_I_Depressed(subject).MeanOccurrence(7)];
            Coverage_G = [Coverage_G outputStats.BP_I_Depressed(subject).Coverage(7)];

            Group{end + 1} = 'BP_I_Depressed';

            Age = [Age France.BP_I.Depressed{subject}.age];
            Gender{end + 1} = France.BP_I.Depressed{subject}.gender;
            SubjectID = [SubjectID France.BP_I.Depressed{subject}.number];

            GAF = [GAF France.BP_I.Depressed{subject}.GAF];
            MADRS = [MADRS France.BP_I.Depressed{subject}.MADRS];
            YMRS = [YMRS France.BP_I.Depressed{subject}.YMRS];

            AD = [AD France.BP_I.Depressed{subject}.meds.AD];
            AP = [AP France.BP_I.Depressed{subject}.meds.AP];
            MS = [MS France.BP_I.Depressed{subject}.meds.MS];
            ANX = [ANX France.BP_I.Depressed{subject}.meds.ANX];
            OTHER = [OTHER France.BP_I.Depressed{subject}.meds.OTHER];
            poly_count = [poly_count France.BP_I.Depressed{subject}.meds.poly_count];
            poly_ge2 = [poly_ge2 France.BP_I.Depressed{subject}.meds.poly_ge2];
            drug_count = [drug_count France.BP_I.Depressed{subject}.meds.drug_count];
            meds_missing = [meds_missing France.BP_I.Depressed{subject}.meds.meds_missing];

        end

        for subject = 1: size(outputStats.BP_I_Euthymic, 2) % then BP_I_Euthymic

            Duration_A = [Duration_A outputStats.BP_I_Euthymic(subject).MeanDuration(1)];
            Occurrence_A = [Occurrence_A outputStats.BP_I_Euthymic(subject).MeanOccurrence(1)];
            Coverage_A = [Coverage_A outputStats.BP_I_Euthymic(subject).Coverage(1)];

            Duration_B = [Duration_B outputStats.BP_I_Euthymic(subject).MeanDuration(2)];
            Occurrence_B = [Occurrence_B outputStats.BP_I_Euthymic(subject).MeanOccurrence(2)];
            Coverage_B = [Coverage_B outputStats.BP_I_Euthymic(subject).Coverage(2)];

            Duration_C = [Duration_C outputStats.BP_I_Euthymic(subject).MeanDuration(3)];
            Occurrence_C = [Occurrence_C outputStats.BP_I_Euthymic(subject).MeanOccurrence(3)];
            Coverage_C = [Coverage_C outputStats.BP_I_Euthymic(subject).Coverage(3)];

            Duration_D = [Duration_D outputStats.BP_I_Euthymic(subject).MeanDuration(4)];
            Occurrence_D = [Occurrence_D outputStats.BP_I_Euthymic(subject).MeanOccurrence(4)];
            Coverage_D = [Coverage_D outputStats.BP_I_Euthymic(subject).Coverage(4)];

            Duration_E = [Duration_E outputStats.BP_I_Euthymic(subject).MeanDuration(5)];
            Occurrence_E = [Occurrence_E outputStats.BP_I_Euthymic(subject).MeanOccurrence(5)];
            Coverage_E = [Coverage_E outputStats.BP_I_Euthymic(subject).Coverage(5)];

            Duration_F = [Duration_F outputStats.BP_I_Euthymic(subject).MeanDuration(6)];
            Occurrence_F = [Occurrence_F outputStats.BP_I_Euthymic(subject).MeanOccurrence(6)];
            Coverage_F = [Coverage_F outputStats.BP_I_Euthymic(subject).Coverage(6)];

            Duration_G = [Duration_G outputStats.BP_I_Euthymic(subject).MeanDuration(7)];
            Occurrence_G = [Occurrence_G outputStats.BP_I_Euthymic(subject).MeanOccurrence(7)];
            Coverage_G = [Coverage_G outputStats.BP_I_Euthymic(subject).Coverage(7)];

            Group{end + 1}  = 'BP_I_Euthymic';

            Age = [Age France.BP_I.Euthymic{subject}.age];
            Gender{end + 1} = France.BP_I.Euthymic{subject}.gender;
            SubjectID = [SubjectID France.BP_I.Euthymic{subject}.number];

            GAF = [GAF France.BP_I.Euthymic{subject}.GAF];
            MADRS = [MADRS France.BP_I.Euthymic{subject}.MADRS];
            YMRS = [YMRS France.BP_I.Euthymic{subject}.YMRS];

            AD = [AD France.BP_I.Euthymic{subject}.meds.AD];
            AP = [AP France.BP_I.Euthymic{subject}.meds.AP];
            MS = [MS France.BP_I.Euthymic{subject}.meds.MS];
            ANX = [ANX France.BP_I.Euthymic{subject}.meds.ANX];
            OTHER = [OTHER France.BP_I.Euthymic{subject}.meds.OTHER];
            poly_count = [poly_count France.BP_I.Euthymic{subject}.meds.poly_count];
            poly_ge2 = [poly_ge2 France.BP_I.Euthymic{subject}.meds.poly_ge2];
            drug_count = [drug_count France.BP_I.Euthymic{subject}.meds.drug_count];
            meds_missing = [meds_missing France.BP_I.Euthymic{subject}.meds.meds_missing];
        end

        for subject = 1: size(outputStats.BP_II_Depressed, 2) % then BP_II_Depressed

            Duration_A = [Duration_A outputStats.BP_II_Depressed(subject).MeanDuration(1)];
            Occurrence_A = [Occurrence_A outputStats.BP_II_Depressed(subject).MeanOccurrence(1)];
            Coverage_A = [Coverage_A outputStats.BP_II_Depressed(subject).Coverage(1)];

            Duration_B = [Duration_B outputStats.BP_II_Depressed(subject).MeanDuration(2)];
            Occurrence_B = [Occurrence_B outputStats.BP_II_Depressed(subject).MeanOccurrence(2)];
            Coverage_B = [Coverage_B outputStats.BP_II_Depressed(subject).Coverage(2)];

            Duration_C = [Duration_C outputStats.BP_II_Depressed(subject).MeanDuration(3)];
            Occurrence_C = [Occurrence_C outputStats.BP_II_Depressed(subject).MeanOccurrence(3)];
            Coverage_C = [Coverage_C outputStats.BP_II_Depressed(subject).Coverage(3)];

            Duration_D = [Duration_D outputStats.BP_II_Depressed(subject).MeanDuration(4)];
            Occurrence_D = [Occurrence_D outputStats.BP_II_Depressed(subject).MeanOccurrence(4)];
            Coverage_D = [Coverage_D outputStats.BP_II_Depressed(subject).Coverage(4)];

            Duration_E = [Duration_E outputStats.BP_II_Depressed(subject).MeanDuration(5)];
            Occurrence_E = [Occurrence_E outputStats.BP_II_Depressed(subject).MeanOccurrence(5)];
            Coverage_E = [Coverage_E outputStats.BP_II_Depressed(subject).Coverage(5)];

            Duration_F = [Duration_F outputStats.BP_II_Depressed(subject).MeanDuration(6)];
            Occurrence_F = [Occurrence_F outputStats.BP_II_Depressed(subject).MeanOccurrence(6)];
            Coverage_F = [Coverage_F outputStats.BP_II_Depressed(subject).Coverage(6)];

            Duration_G = [Duration_G outputStats.BP_II_Depressed(subject).MeanDuration(7)];
            Occurrence_G = [Occurrence_G outputStats.BP_II_Depressed(subject).MeanOccurrence(7)];
            Coverage_G = [Coverage_G outputStats.BP_II_Depressed(subject).Coverage(7)];

            Group{end + 1} = 'BP_II_Depressed';

            Age = [Age France.BP_II.Depressed{subject}.age];
            Gender{end + 1} = France.BP_II.Depressed{subject}.gender;
            SubjectID = [SubjectID France.BP_II.Depressed{subject}.number];

            GAF = [GAF France.BP_II.Depressed{subject}.GAF];
            MADRS = [MADRS France.BP_II.Depressed{subject}.MADRS];
            YMRS = [YMRS France.BP_II.Depressed{subject}.YMRS];

            AD = [AD France.BP_II.Depressed{subject}.meds.AD];
            AP = [AP France.BP_II.Depressed{subject}.meds.AP];
            MS = [MS France.BP_II.Depressed{subject}.meds.MS];
            ANX = [ANX France.BP_II.Depressed{subject}.meds.ANX];
            OTHER = [OTHER France.BP_II.Depressed{subject}.meds.OTHER];
            poly_count = [poly_count France.BP_II.Depressed{subject}.meds.poly_count];
            poly_ge2 = [poly_ge2 France.BP_II.Depressed{subject}.meds.poly_ge2];
            drug_count = [drug_count France.BP_II.Depressed{subject}.meds.drug_count];
            meds_missing = [meds_missing France.BP_II.Depressed{subject}.meds.meds_missing];
        end

        for subject = 1: size(outputStats.BP_II_Euthymic, 2) % then BP_II_Euthymic

            Duration_A = [Duration_A outputStats.BP_II_Euthymic(subject).MeanDuration(1)];
            Occurrence_A = [Occurrence_A outputStats.BP_II_Euthymic(subject).MeanOccurrence(1)];
            Coverage_A = [Coverage_A outputStats.BP_II_Euthymic(subject).Coverage(1)];

            Duration_B = [Duration_B outputStats.BP_II_Euthymic(subject).MeanDuration(2)];
            Occurrence_B = [Occurrence_B outputStats.BP_II_Euthymic(subject).MeanOccurrence(2)];
            Coverage_B = [Coverage_B outputStats.BP_II_Euthymic(subject).Coverage(2)];

            Duration_C = [Duration_C outputStats.BP_II_Euthymic(subject).MeanDuration(3)];
            Occurrence_C = [Occurrence_C outputStats.BP_II_Euthymic(subject).MeanOccurrence(3)];
            Coverage_C = [Coverage_C outputStats.BP_II_Euthymic(subject).Coverage(3)];

            Duration_D = [Duration_D outputStats.BP_II_Euthymic(subject).MeanDuration(4)];
            Occurrence_D = [Occurrence_D outputStats.BP_II_Euthymic(subject).MeanOccurrence(4)];
            Coverage_D = [Coverage_D outputStats.BP_II_Euthymic(subject).Coverage(4)];

            Duration_E = [Duration_E outputStats.BP_II_Euthymic(subject).MeanDuration(5)];
            Occurrence_E = [Occurrence_E outputStats.BP_II_Euthymic(subject).MeanOccurrence(5)];
            Coverage_E = [Coverage_E outputStats.BP_II_Euthymic(subject).Coverage(5)];

            Duration_F = [Duration_F outputStats.BP_II_Euthymic(subject).MeanDuration(6)];
            Occurrence_F = [Occurrence_F outputStats.BP_II_Euthymic(subject).MeanOccurrence(6)];
            Coverage_F = [Coverage_F outputStats.BP_II_Euthymic(subject).Coverage(6)];

            Duration_G = [Duration_G outputStats.BP_II_Euthymic(subject).MeanDuration(7)];
            Occurrence_G = [Occurrence_G outputStats.BP_II_Euthymic(subject).MeanOccurrence(7)];
            Coverage_G = [Coverage_G outputStats.BP_II_Euthymic(subject).Coverage(7)];

            Group{end + 1} = 'BP_II_Euthymic';

            Age = [Age France.BP_II.Euthymic{subject}.age];
            Gender{end + 1} = France.BP_II.Euthymic{subject}.gender;
            SubjectID = [SubjectID France.BP_II.Euthymic{subject}.number];

            GAF = [GAF France.BP_II.Euthymic{subject}.GAF];
            MADRS = [MADRS France.BP_II.Euthymic{subject}.MADRS];
            YMRS = [YMRS France.BP_II.Euthymic{subject}.YMRS];

            AD = [AD France.BP_II.Euthymic{subject}.meds.AD];
            AP = [AP France.BP_II.Euthymic{subject}.meds.AP];
            MS = [MS France.BP_II.Euthymic{subject}.meds.MS];
            ANX = [ANX France.BP_II.Euthymic{subject}.meds.ANX];
            OTHER = [OTHER France.BP_II.Euthymic{subject}.meds.OTHER];
            poly_count = [poly_count France.BP_II.Euthymic{subject}.meds.poly_count];
            poly_ge2 = [poly_ge2 France.BP_II.Euthymic{subject}.meds.poly_ge2];
            drug_count = [drug_count France.BP_II.Euthymic{subject}.meds.drug_count];
            meds_missing = [meds_missing France.BP_II.Euthymic{subject}.meds.meds_missing];

        end

        for subject = 1: size(outputStats.HC, 2) % then HC

            Duration_A = [Duration_A outputStats.HC(subject).MeanDuration(1)];
            Occurrence_A = [Occurrence_A outputStats.HC(subject).MeanOccurrence(1)];
            Coverage_A = [Coverage_A outputStats.HC(subject).Coverage(1)];

            Duration_B = [Duration_B outputStats.HC(subject).MeanDuration(2)];
            Occurrence_B = [Occurrence_B outputStats.HC(subject).MeanOccurrence(2)];
            Coverage_B = [Coverage_B outputStats.HC(subject).Coverage(2)];

            Duration_C = [Duration_C outputStats.HC(subject).MeanDuration(3)];
            Occurrence_C = [Occurrence_C outputStats.HC(subject).MeanOccurrence(3)];
            Coverage_C = [Coverage_C outputStats.HC(subject).Coverage(3)];

            Duration_D = [Duration_D outputStats.HC(subject).MeanDuration(4)];
            Occurrence_D = [Occurrence_D outputStats.HC(subject).MeanOccurrence(4)];
            Coverage_D = [Coverage_D outputStats.HC(subject).Coverage(4)];

            Duration_E = [Duration_E outputStats.HC(subject).MeanDuration(5)];
            Occurrence_E = [Occurrence_E outputStats.HC(subject).MeanOccurrence(5)];
            Coverage_E = [Coverage_E outputStats.HC(subject).Coverage(5)];

            Duration_F = [Duration_F outputStats.HC(subject).MeanDuration(6)];
            Occurrence_F = [Occurrence_F outputStats.HC(subject).MeanOccurrence(6)];
            Coverage_F = [Coverage_F outputStats.HC(subject).Coverage(6)];

            Duration_G = [Duration_G outputStats.HC(subject).MeanDuration(7)];
            Occurrence_G = [Occurrence_G outputStats.HC(subject).MeanOccurrence(7)];
            Coverage_G = [Coverage_G outputStats.HC(subject).Coverage(7)];

            Group{end + 1}  = 'HC';

            Age = [Age France.HC{subject}.age];
            Gender{end + 1} = France.HC{subject}.gender;
            SubjectID = [SubjectID France.HC{subject}.number];

            GAF = [GAF France.HC{subject}.GAF];
            MADRS = [MADRS France.HC{subject}.MADRS];
            YMRS = [YMRS France.HC{subject}.YMRS];

            AD = [AD 0];
            AP = [AP 0];
            MS = [MS 0];
            ANX = [ANX 0];
            OTHER = [OTHER 0];
            poly_count = [poly_count 0];
            poly_ge2 = [poly_ge2 0];
            drug_count = [drug_count 0];
            meds_missing = [meds_missing 0];

        end

        for subject = 1: size(outputStats.Siblings, 2) % then siblings

            Duration_A = [Duration_A outputStats.Siblings(subject).MeanDuration(1)];
            Occurrence_A = [Occurrence_A outputStats.Siblings(subject).MeanOccurrence(1)];
            Coverage_A = [Coverage_A outputStats.Siblings(subject).Coverage(1)];

            Duration_B = [Duration_B outputStats.Siblings(subject).MeanDuration(2)];
            Occurrence_B = [Occurrence_B outputStats.Siblings(subject).MeanOccurrence(2)];
            Coverage_B = [Coverage_B outputStats.Siblings(subject).Coverage(2)];

            Duration_C = [Duration_C outputStats.Siblings(subject).MeanDuration(3)];
            Occurrence_C = [Occurrence_C outputStats.Siblings(subject).MeanOccurrence(3)];
            Coverage_C = [Coverage_C outputStats.Siblings(subject).Coverage(3)];

            Duration_D = [Duration_D outputStats.Siblings(subject).MeanDuration(4)];
            Occurrence_D = [Occurrence_D outputStats.Siblings(subject).MeanOccurrence(4)];
            Coverage_D = [Coverage_D outputStats.Siblings(subject).Coverage(4)];

            Duration_E = [Duration_E outputStats.Siblings(subject).MeanDuration(5)];
            Occurrence_E = [Occurrence_E outputStats.Siblings(subject).MeanOccurrence(5)];
            Coverage_E = [Coverage_E outputStats.Siblings(subject).Coverage(5)];

            Duration_F = [Duration_F outputStats.Siblings(subject).MeanDuration(6)];
            Occurrence_F = [Occurrence_F outputStats.Siblings(subject).MeanOccurrence(6)];
            Coverage_F = [Coverage_F outputStats.Siblings(subject).Coverage(6)];

            Duration_G = [Duration_G outputStats.Siblings(subject).MeanDuration(7)];
            Occurrence_G = [Occurrence_G outputStats.Siblings(subject).MeanOccurrence(7)];
            Coverage_G = [Coverage_G outputStats.Siblings(subject).Coverage(7)];

            Group{end + 1}  = 'Siblings';

            Age = [Age France.Siblings{subject}.age];
            Gender{end + 1} = France.Siblings{subject}.gender;
            SubjectID = [SubjectID France.Siblings{subject}.number];

            GAF = [GAF France.Siblings{subject}.GAF];
            MADRS = [MADRS France.Siblings{subject}.MADRS];
            YMRS = [YMRS France.Siblings{subject}.YMRS];

            AD = [AD 0];
            AP = [AP 0];
            MS = [MS 0];
            ANX = [ANX 0];
            OTHER = [OTHER 0];
            poly_count = [poly_count 0];
            poly_ge2 = [poly_ge2 0];
            drug_count = [drug_count 0];
            meds_missing = [meds_missing 0];

        end

        Class_A_BP_I_Depressed_Duration = [];       Class_B_BP_I_Depressed_Duration = [];        Class_C_BP_I_Depressed_Duration = [];          Class_D_BP_I_Depressed_Duration = [];          Class_E_BP_I_Depressed_Duration = [];            Class_F_BP_I_Depressed_Duration = [];        Class_G_BP_I_Depressed_Duration = [];
        Class_A_BP_I_Euthymic_Duration = [];          Class_B_BP_I_Euthymic_Duration = [];          Class_C_BP_I_Euthymic_Duration = [];             Class_D_BP_I_Euthymic_Duration = [];            Class_E_BP_I_Euthymic_Duration = [];               Class_F_BP_I_Euthymic_Duration = [];           Class_G_BP_I_Euthymic_Duration = [];
        Class_A_BP_II_Depressed_Duration = [];      Class_B_BP_II_Depressed_Duration = [];       Class_C_BP_II_Depressed_Duration = [];         Class_D_BP_II_Depressed_Duration = [];         Class_E_BP_II_Depressed_Duration = [];           Class_F_BP_II_Depressed_Duration = [];       Class_G_BP_II_Depressed_Duration = [];
        Class_A_BP_II_Euthymic_Duration = [];         Class_B_BP_II_Euthymic_Duration = [];         Class_C_BP_II_Euthymic_Duration = [];            Class_D_BP_II_Euthymic_Duration = [];            Class_E_BP_II_Euthymic_Duration = [];             Class_F_BP_II_Euthymic_Duration = [];          Class_G_BP_II_Euthymic_Duration = [];
        Class_A_HC_Duration = [];                             Class_B_HC_Duration = [];                              Class_C_HC_Duration = [];                                Class_D_HC_Duration = [];                                Class_E_HC_Duration = [];                                  Class_F_HC_Duration = [];                              Class_G_HC_Duration = [];
        Class_A_Siblings_Duration = [];                      Class_B_Siblings_Duration = [];                       Class_C_Siblings_Duration = [];                        Class_D_Siblings_Duration= [];                          Class_E_Siblings_Duration = [];                          Class_F_Siblings_Duration = [];                       Class_G_Siblings_Duration = [];

        Class_A_BP_I_Depressed_Occurrence = [];  Class_B_BP_I_Depressed_Occurrence = [];    Class_C_BP_I_Depressed_Occurrence = [];    Class_D_BP_I_Depressed_Occurrence = [];     Class_E_BP_I_Depressed_Occurrence = [];      Class_F_BP_I_Depressed_Occurrence = [];    Class_G_BP_I_Depressed_Occurrence = [];
        Class_A_BP_I_Euthymic_Occurrence = [];     Class_B_BP_I_Euthymic_Occurrence = [];      Class_C_BP_I_Euthymic_Occurrence = [];       Class_D_BP_I_Euthymic_Occurrence = [];        Class_E_BP_I_Euthymic_Occurrence = [];        Class_F_BP_I_Euthymic_Occurrence = [];       Class_G_BP_I_Euthymic_Occurrence = [];
        Class_A_BP_II_Depressed_Occurrence = []; Class_B_BP_II_Depressed_Occurrence = [];   Class_C_BP_II_Depressed_Occurrence = [];   Class_D_BP_II_Depressed_Occurrence = [];     Class_E_BP_II_Depressed_Occurrence = [];    Class_F_BP_II_Depressed_Occurrence = [];   Class_G_BP_II_Depressed_Occurrence = [];
        Class_A_BP_II_Euthymic_Occurrence = [];    Class_B_BP_II_Euthymic_Occurrence = [];     Class_C_BP_II_Euthymic_Occurrence = [];      Class_D_BP_II_Euthymic_Occurrence = [];       Class_E_BP_II_Euthymic_Occurrence = [];       Class_F_BP_II_Euthymic_Occurrence = [];      Class_G_BP_II_Euthymic_Occurrence = [];
        Class_A_HC_Occurrence = [];                        Class_B_HC_Occurrence = [];                          Class_C_HC_Occurrence = [];                          Class_D_HC_Occurrence = [];                            Class_E_HC_Occurrence = [];                           Class_F_HC_Occurrence = [];                          Class_G_HC_Occurrence = [];
        Class_A_Siblings_Occurrence = [];                 Class_B_Siblings_Occurrence = [];                   Class_C_Siblings_Occurrence = [];                  Class_D_Siblings_Occurrence= [];                     Class_E_Siblings_Occurrence = [];                    Class_F_Siblings_Occurrence = [];                   Class_G_Siblings_Occurrence = [];

        Class_A_BP_I_Depressed_Coverage = [];      Class_B_BP_I_Depressed_Coverage = [];      Class_C_BP_I_Depressed_Coverage = [];       Class_D_BP_I_Depressed_Coverage = [];        Class_E_BP_I_Depressed_Coverage = [];        Class_F_BP_I_Depressed_Coverage = [];       Class_G_BP_I_Depressed_Coverage = [];
        Class_A_BP_I_Euthymic_Coverage = [];         Class_B_BP_I_Euthymic_Coverage = [];        Class_C_BP_I_Euthymic_Coverage = [];          Class_D_BP_I_Euthymic_Coverage = [];          Class_E_BP_I_Euthymic_Coverage = [];           Class_F_BP_I_Euthymic_Coverage = [];          Class_G_BP_I_Euthymic_Coverage = [];
        Class_A_BP_II_Depressed_Coverage = [];      Class_B_BP_II_Depressed_Coverage = [];    Class_C_BP_II_Depressed_Coverage = [];      Class_D_BP_II_Depressed_Coverage = [];       Class_E_BP_II_Depressed_Coverage = [];       Class_F_BP_II_Depressed_Coverage = [];      Class_G_BP_II_Depressed_Coverage = [];
        Class_A_BP_II_Euthymic_Coverage = [];         Class_B_BP_II_Euthymic_Coverage = [];      Class_C_BP_II_Euthymic_Coverage = [];         Class_D_BP_II_Euthymic_Coverage = [];         Class_E_BP_II_Euthymic_Coverage = [];          Class_F_BP_II_Euthymic_Coverage = [];         Class_G_BP_II_Euthymic_Coverage = [];
        Class_A_HC_Coverage = [];                             Class_B_HC_Coverage = [];                           Class_C_HC_Coverage = [];                             Class_D_HC_Coverage = [];                              Class_E_HC_Coverage = [];                              Class_F_HC_Coverage = [];                             Class_G_HC_Coverage = [];
        Class_A_Siblings_Coverage = [];                      Class_B_Siblings_Coverage = [];                   Class_C_Siblings_Coverage = [];                      Class_D_Siblings_Coverage= [];                       Class_E_Siblings_Coverage = [];                       Class_F_Siblings_Coverage = [];                      Class_G_Siblings_Coverage = [];

        % Now, pool the three metrics of each microstate across all subjects in a group

        for subject = 1: size(outputStats.BP_I_Depressed, 2) % First, BP_I_Depressed

            Class_A_BP_I_Depressed_Duration      = [Class_A_BP_I_Depressed_Duration      outputStats.BP_I_Depressed(subject).MeanDuration(1)];
            Class_A_BP_I_Depressed_Occurrence = [Class_A_BP_I_Depressed_Occurrence outputStats.BP_I_Depressed(subject).MeanOccurrence(1)];
            Class_A_BP_I_Depressed_Coverage    = [Class_A_BP_I_Depressed_Coverage    outputStats.BP_I_Depressed(subject).Coverage(1)];

            Class_B_BP_I_Depressed_Duration      = [Class_B_BP_I_Depressed_Duration      outputStats.BP_I_Depressed(subject).MeanDuration(2)];
            Class_B_BP_I_Depressed_Occurrence = [Class_B_BP_I_Depressed_Occurrence outputStats.BP_I_Depressed(subject).MeanOccurrence(2)];
            Class_B_BP_I_Depressed_Coverage    = [Class_B_BP_I_Depressed_Coverage    outputStats.BP_I_Depressed(subject).Coverage(2)];

            Class_C_BP_I_Depressed_Duration      = [Class_C_BP_I_Depressed_Duration      outputStats.BP_I_Depressed(subject).MeanDuration(3)];
            Class_C_BP_I_Depressed_Occurrence = [Class_C_BP_I_Depressed_Occurrence outputStats.BP_I_Depressed(subject).MeanOccurrence(3)];
            Class_C_BP_I_Depressed_Coverage    = [Class_C_BP_I_Depressed_Coverage    outputStats.BP_I_Depressed(subject).Coverage(3)];

            Class_D_BP_I_Depressed_Duration      = [Class_D_BP_I_Depressed_Duration      outputStats.BP_I_Depressed(subject).MeanDuration(4)];
            Class_D_BP_I_Depressed_Occurrence = [Class_D_BP_I_Depressed_Occurrence outputStats.BP_I_Depressed(subject).MeanOccurrence(4)];
            Class_D_BP_I_Depressed_Coverage    = [Class_D_BP_I_Depressed_Coverage    outputStats.BP_I_Depressed(subject).Coverage(4)];

            Class_E_BP_I_Depressed_Duration      = [Class_E_BP_I_Depressed_Duration      outputStats.BP_I_Depressed(subject).MeanDuration(5)];
            Class_E_BP_I_Depressed_Occurrence = [Class_E_BP_I_Depressed_Occurrence outputStats.BP_I_Depressed(subject).MeanOccurrence(5)];
            Class_E_BP_I_Depressed_Coverage    = [Class_E_BP_I_Depressed_Coverage    outputStats.BP_I_Depressed(subject).Coverage(5)];

            Class_F_BP_I_Depressed_Duration      = [Class_F_BP_I_Depressed_Duration      outputStats.BP_I_Depressed(subject).MeanDuration(6)];
            Class_F_BP_I_Depressed_Occurrence = [Class_F_BP_I_Depressed_Occurrence outputStats.BP_I_Depressed(subject).MeanOccurrence(6)];
            Class_F_BP_I_Depressed_Coverage    = [Class_F_BP_I_Depressed_Coverage    outputStats.BP_I_Depressed(subject).Coverage(6)];

            Class_G_BP_I_Depressed_Duration      = [Class_G_BP_I_Depressed_Duration      outputStats.BP_I_Depressed(subject).MeanDuration(7)];
            Class_G_BP_I_Depressed_Occurrence = [Class_G_BP_I_Depressed_Occurrence outputStats.BP_I_Depressed(subject).MeanOccurrence(7)];
            Class_G_BP_I_Depressed_Coverage    = [Class_G_BP_I_Depressed_Coverage    outputStats.BP_I_Depressed(subject).Coverage(7)];

        end

        for subject = 1: size(outputStats.BP_I_Euthymic, 2) % then, BP_I_Euthymic

            Class_A_BP_I_Euthymic_Duration      = [Class_A_BP_I_Euthymic_Duration      outputStats.BP_I_Euthymic(subject).MeanDuration(1)];
            Class_A_BP_I_Euthymic_Occurrence = [Class_A_BP_I_Euthymic_Occurrence outputStats.BP_I_Euthymic(subject).MeanOccurrence(1)];
            Class_A_BP_I_Euthymic_Coverage    = [Class_A_BP_I_Euthymic_Coverage    outputStats.BP_I_Euthymic(subject).Coverage(1)];

            Class_B_BP_I_Euthymic_Duration      = [Class_B_BP_I_Euthymic_Duration      outputStats.BP_I_Euthymic(subject).MeanDuration(2)];
            Class_B_BP_I_Euthymic_Occurrence = [Class_B_BP_I_Euthymic_Occurrence outputStats.BP_I_Euthymic(subject).MeanOccurrence(2)];
            Class_B_BP_I_Euthymic_Coverage    = [Class_B_BP_I_Euthymic_Coverage    outputStats.BP_I_Euthymic(subject).Coverage(2)];

            Class_C_BP_I_Euthymic_Duration      = [Class_C_BP_I_Euthymic_Duration      outputStats.BP_I_Euthymic(subject).MeanDuration(3)];
            Class_C_BP_I_Euthymic_Occurrence = [Class_C_BP_I_Euthymic_Occurrence outputStats.BP_I_Euthymic(subject).MeanOccurrence(3)];
            Class_C_BP_I_Euthymic_Coverage    = [Class_C_BP_I_Euthymic_Coverage    outputStats.BP_I_Euthymic(subject).Coverage(3)];

            Class_D_BP_I_Euthymic_Duration      = [Class_D_BP_I_Euthymic_Duration      outputStats.BP_I_Euthymic(subject).MeanDuration(4)];
            Class_D_BP_I_Euthymic_Occurrence = [Class_D_BP_I_Euthymic_Occurrence outputStats.BP_I_Euthymic(subject).MeanOccurrence(4)];
            Class_D_BP_I_Euthymic_Coverage    = [Class_D_BP_I_Euthymic_Coverage    outputStats.BP_I_Euthymic(subject).Coverage(4)];

            Class_E_BP_I_Euthymic_Duration      = [Class_E_BP_I_Euthymic_Duration      outputStats.BP_I_Euthymic(subject).MeanDuration(5)];
            Class_E_BP_I_Euthymic_Occurrence = [Class_E_BP_I_Euthymic_Occurrence outputStats.BP_I_Euthymic(subject).MeanOccurrence(5)];
            Class_E_BP_I_Euthymic_Coverage    = [Class_E_BP_I_Euthymic_Coverage    outputStats.BP_I_Euthymic(subject).Coverage(5)];

            Class_F_BP_I_Euthymic_Duration      = [Class_F_BP_I_Euthymic_Duration      outputStats.BP_I_Euthymic(subject).MeanDuration(6)];
            Class_F_BP_I_Euthymic_Occurrence = [Class_F_BP_I_Euthymic_Occurrence outputStats.BP_I_Euthymic(subject).MeanOccurrence(6)];
            Class_F_BP_I_Euthymic_Coverage    = [Class_F_BP_I_Euthymic_Coverage    outputStats.BP_I_Euthymic(subject).Coverage(6)];

            Class_G_BP_I_Euthymic_Duration      = [Class_G_BP_I_Euthymic_Duration      outputStats.BP_I_Euthymic(subject).MeanDuration(7)];
            Class_G_BP_I_Euthymic_Occurrence = [Class_G_BP_I_Euthymic_Occurrence outputStats.BP_I_Euthymic(subject).MeanOccurrence(7)];
            Class_G_BP_I_Euthymic_Coverage    = [Class_G_BP_I_Euthymic_Coverage    outputStats.BP_I_Euthymic(subject).Coverage(7)];

        end

        for subject = 1: size(outputStats.BP_II_Depressed, 2) % then, BP_II_Depressed

            Class_A_BP_II_Depressed_Duration      = [Class_A_BP_II_Depressed_Duration      outputStats.BP_II_Depressed(subject).MeanDuration(1)];
            Class_A_BP_II_Depressed_Occurrence = [Class_A_BP_II_Depressed_Occurrence outputStats.BP_II_Depressed(subject).MeanOccurrence(1)];
            Class_A_BP_II_Depressed_Coverage    = [Class_A_BP_II_Depressed_Coverage    outputStats.BP_II_Depressed(subject).Coverage(1)];

            Class_B_BP_II_Depressed_Duration      = [Class_B_BP_II_Depressed_Duration      outputStats.BP_II_Depressed(subject).MeanDuration(2)];
            Class_B_BP_II_Depressed_Occurrence = [Class_B_BP_II_Depressed_Occurrence outputStats.BP_II_Depressed(subject).MeanOccurrence(2)];
            Class_B_BP_II_Depressed_Coverage    = [Class_B_BP_II_Depressed_Coverage    outputStats.BP_II_Depressed(subject).Coverage(2)];

            Class_C_BP_II_Depressed_Duration      = [Class_C_BP_II_Depressed_Duration      outputStats.BP_II_Depressed(subject).MeanDuration(3)];
            Class_C_BP_II_Depressed_Occurrence = [Class_C_BP_II_Depressed_Occurrence outputStats.BP_II_Depressed(subject).MeanOccurrence(3)];
            Class_C_BP_II_Depressed_Coverage    = [Class_C_BP_II_Depressed_Coverage    outputStats.BP_II_Depressed(subject).Coverage(3)];

            Class_D_BP_II_Depressed_Duration      = [Class_D_BP_II_Depressed_Duration      outputStats.BP_II_Depressed(subject).MeanDuration(4)];
            Class_D_BP_II_Depressed_Occurrence = [Class_D_BP_II_Depressed_Occurrence outputStats.BP_II_Depressed(subject).MeanOccurrence(4)];
            Class_D_BP_II_Depressed_Coverage    = [Class_D_BP_II_Depressed_Coverage    outputStats.BP_II_Depressed(subject).Coverage(4)];

            Class_E_BP_II_Depressed_Duration      = [Class_E_BP_II_Depressed_Duration      outputStats.BP_II_Depressed(subject).MeanDuration(5)];
            Class_E_BP_II_Depressed_Occurrence = [Class_E_BP_II_Depressed_Occurrence outputStats.BP_II_Depressed(subject).MeanOccurrence(5)];
            Class_E_BP_II_Depressed_Coverage    = [Class_E_BP_II_Depressed_Coverage    outputStats.BP_II_Depressed(subject).Coverage(5)];

            Class_F_BP_II_Depressed_Duration      = [Class_F_BP_II_Depressed_Duration      outputStats.BP_II_Depressed(subject).MeanDuration(6)];
            Class_F_BP_II_Depressed_Occurrence = [Class_F_BP_II_Depressed_Occurrence outputStats.BP_II_Depressed(subject).MeanOccurrence(6)];
            Class_F_BP_II_Depressed_Coverage    = [Class_F_BP_II_Depressed_Coverage    outputStats.BP_II_Depressed(subject).Coverage(6)];

            Class_G_BP_II_Depressed_Duration      = [Class_G_BP_II_Depressed_Duration      outputStats.BP_II_Depressed(subject).MeanDuration(7)];
            Class_G_BP_II_Depressed_Occurrence = [Class_G_BP_II_Depressed_Occurrence outputStats.BP_II_Depressed(subject).MeanOccurrence(7)];
            Class_G_BP_II_Depressed_Coverage    = [Class_G_BP_II_Depressed_Coverage    outputStats.BP_II_Depressed(subject).Coverage(7)];

        end

        for subject = 1: size(outputStats.BP_II_Euthymic, 2) % then, BP_II_Euthymic

            Class_A_BP_II_Euthymic_Duration      = [Class_A_BP_II_Euthymic_Duration      outputStats.BP_II_Euthymic(subject).MeanDuration(1)];
            Class_A_BP_II_Euthymic_Occurrence = [Class_A_BP_II_Euthymic_Occurrence outputStats.BP_II_Euthymic(subject).MeanOccurrence(1)];
            Class_A_BP_II_Euthymic_Coverage    = [Class_A_BP_II_Euthymic_Coverage    outputStats.BP_II_Euthymic(subject).Coverage(1)];

            Class_B_BP_II_Euthymic_Duration      = [Class_B_BP_II_Euthymic_Duration      outputStats.BP_II_Euthymic(subject).MeanDuration(2)];
            Class_B_BP_II_Euthymic_Occurrence = [Class_B_BP_II_Euthymic_Occurrence outputStats.BP_II_Euthymic(subject).MeanOccurrence(2)];
            Class_B_BP_II_Euthymic_Coverage    = [Class_B_BP_II_Euthymic_Coverage    outputStats.BP_II_Euthymic(subject).Coverage(2)];

            Class_C_BP_II_Euthymic_Duration      = [Class_C_BP_II_Euthymic_Duration      outputStats.BP_II_Euthymic(subject).MeanDuration(3)];
            Class_C_BP_II_Euthymic_Occurrence = [Class_C_BP_II_Euthymic_Occurrence outputStats.BP_II_Euthymic(subject).MeanOccurrence(3)];
            Class_C_BP_II_Euthymic_Coverage    = [Class_C_BP_II_Euthymic_Coverage    outputStats.BP_II_Euthymic(subject).Coverage(3)];

            Class_D_BP_II_Euthymic_Duration      = [Class_D_BP_II_Euthymic_Duration      outputStats.BP_II_Euthymic(subject).MeanDuration(4)];
            Class_D_BP_II_Euthymic_Occurrence = [Class_D_BP_II_Euthymic_Occurrence outputStats.BP_II_Euthymic(subject).MeanOccurrence(4)];
            Class_D_BP_II_Euthymic_Coverage    = [Class_D_BP_II_Euthymic_Coverage    outputStats.BP_II_Euthymic(subject).Coverage(4)];

            Class_E_BP_II_Euthymic_Duration      = [Class_E_BP_II_Euthymic_Duration      outputStats.BP_II_Euthymic(subject).MeanDuration(5)];
            Class_E_BP_II_Euthymic_Occurrence = [Class_E_BP_II_Euthymic_Occurrence outputStats.BP_II_Euthymic(subject).MeanOccurrence(5)];
            Class_E_BP_II_Euthymic_Coverage    = [Class_E_BP_II_Euthymic_Coverage    outputStats.BP_II_Euthymic(subject).Coverage(5)];

            Class_F_BP_II_Euthymic_Duration      = [Class_F_BP_II_Euthymic_Duration      outputStats.BP_II_Euthymic(subject).MeanDuration(6)];
            Class_F_BP_II_Euthymic_Occurrence = [Class_F_BP_II_Euthymic_Occurrence outputStats.BP_II_Euthymic(subject).MeanOccurrence(6)];
            Class_F_BP_II_Euthymic_Coverage    = [Class_F_BP_II_Euthymic_Coverage    outputStats.BP_II_Euthymic(subject).Coverage(6)];

            Class_G_BP_II_Euthymic_Duration      = [Class_G_BP_II_Euthymic_Duration      outputStats.BP_II_Euthymic(subject).MeanDuration(7)];
            Class_G_BP_II_Euthymic_Occurrence = [Class_G_BP_II_Euthymic_Occurrence outputStats.BP_II_Euthymic(subject).MeanOccurrence(7)];
            Class_G_BP_II_Euthymic_Coverage    = [Class_G_BP_II_Euthymic_Coverage    outputStats.BP_II_Euthymic(subject).Coverage(7)];

        end

        for subject = 1: size(outputStats.HC, 2) % then, HC

            Class_A_HC_Duration      = [Class_A_HC_Duration      outputStats.HC(subject).MeanDuration(1)];
            Class_A_HC_Occurrence = [Class_A_HC_Occurrence outputStats.HC(subject).MeanOccurrence(1)];
            Class_A_HC_Coverage    = [Class_A_HC_Coverage    outputStats.HC(subject).Coverage(1)];

            Class_B_HC_Duration      = [Class_B_HC_Duration      outputStats.HC(subject).MeanDuration(2)];
            Class_B_HC_Occurrence = [Class_B_HC_Occurrence outputStats.HC(subject).MeanOccurrence(2)];
            Class_B_HC_Coverage    = [Class_B_HC_Coverage    outputStats.HC(subject).Coverage(2)];

            Class_C_HC_Duration      = [Class_C_HC_Duration      outputStats.HC(subject).MeanDuration(3)];
            Class_C_HC_Occurrence = [Class_C_HC_Occurrence outputStats.HC(subject).MeanOccurrence(3)];
            Class_C_HC_Coverage    = [Class_C_HC_Coverage    outputStats.HC(subject).Coverage(3)];

            Class_D_HC_Duration      = [Class_D_HC_Duration      outputStats.HC(subject).MeanDuration(4)];
            Class_D_HC_Occurrence = [Class_D_HC_Occurrence outputStats.HC(subject).MeanOccurrence(4)];
            Class_D_HC_Coverage    = [Class_D_HC_Coverage    outputStats.HC(subject).Coverage(4)];

            Class_E_HC_Duration      = [Class_E_HC_Duration      outputStats.HC(subject).MeanDuration(5)];
            Class_E_HC_Occurrence = [Class_E_HC_Occurrence outputStats.HC(subject).MeanOccurrence(5)];
            Class_E_HC_Coverage    = [Class_E_HC_Coverage    outputStats.HC(subject).Coverage(5)];

            Class_F_HC_Duration      = [Class_F_HC_Duration      outputStats.HC(subject).MeanDuration(6)];
            Class_F_HC_Occurrence = [Class_F_HC_Occurrence outputStats.HC(subject).MeanOccurrence(6)];
            Class_F_HC_Coverage    = [Class_F_HC_Coverage    outputStats.HC(subject).Coverage(6)];

            Class_G_HC_Duration      = [Class_G_HC_Duration      outputStats.HC(subject).MeanDuration(7)];
            Class_G_HC_Occurrence = [Class_G_HC_Occurrence outputStats.HC(subject).MeanOccurrence(7)];
            Class_G_HC_Coverage    = [Class_G_HC_Coverage    outputStats.HC(subject).Coverage(7)];

        end

        for subject = 1: size(outputStats.Siblings, 2) % then, Siblings

            Class_A_Siblings_Duration      = [Class_A_Siblings_Duration      outputStats.Siblings(subject).MeanDuration(1)];
            Class_A_Siblings_Occurrence = [Class_A_Siblings_Occurrence outputStats.Siblings(subject).MeanOccurrence(1)];
            Class_A_Siblings_Coverage    = [Class_A_Siblings_Coverage    outputStats.Siblings(subject).Coverage(1)];

            Class_B_Siblings_Duration      = [Class_B_Siblings_Duration      outputStats.Siblings(subject).MeanDuration(2)];
            Class_B_Siblings_Occurrence = [Class_B_Siblings_Occurrence outputStats.Siblings(subject).MeanOccurrence(2)];
            Class_B_Siblings_Coverage    = [Class_B_Siblings_Coverage    outputStats.Siblings(subject).Coverage(2)];

            Class_C_Siblings_Duration      = [Class_C_Siblings_Duration      outputStats.Siblings(subject).MeanDuration(3)];
            Class_C_Siblings_Occurrence = [Class_C_Siblings_Occurrence outputStats.Siblings(subject).MeanOccurrence(3)];
            Class_C_Siblings_Coverage    = [Class_C_Siblings_Coverage    outputStats.Siblings(subject).Coverage(3)];

            Class_D_Siblings_Duration      = [Class_D_Siblings_Duration      outputStats.Siblings(subject).MeanDuration(4)];
            Class_D_Siblings_Occurrence = [Class_D_Siblings_Occurrence outputStats.Siblings(subject).MeanOccurrence(4)];
            Class_D_Siblings_Coverage    = [Class_D_Siblings_Coverage    outputStats.Siblings(subject).Coverage(4)];

            Class_E_Siblings_Duration      = [Class_E_Siblings_Duration      outputStats.Siblings(subject).MeanDuration(5)];
            Class_E_Siblings_Occurrence = [Class_E_Siblings_Occurrence outputStats.Siblings(subject).MeanOccurrence(5)];
            Class_E_Siblings_Coverage    = [Class_E_Siblings_Coverage    outputStats.Siblings(subject).Coverage(5)];

            Class_F_Siblings_Duration      = [Class_F_Siblings_Duration      outputStats.Siblings(subject).MeanDuration(6)];
            Class_F_Siblings_Occurrence = [Class_F_Siblings_Occurrence outputStats.Siblings(subject).MeanOccurrence(6)];
            Class_F_Siblings_Coverage    = [Class_F_Siblings_Coverage    outputStats.Siblings(subject).Coverage(6)];

            Class_G_Siblings_Duration      = [Class_G_Siblings_Duration      outputStats.Siblings(subject).MeanDuration(7)];
            Class_G_Siblings_Occurrence = [Class_G_Siblings_Occurrence outputStats.Siblings(subject).MeanOccurrence(7)];
            Class_G_Siblings_Coverage    = [Class_G_Siblings_Coverage    outputStats.Siblings(subject).Coverage(7)];

        end

        %% Test for significant effect of group on each of the templates' duration, occurrence and coverage

        clear T

        T = table(SubjectID(:), categorical(Group(:)), Age(:), categorical(Gender(:)), ...
            Duration_A(:) * 1000, Occurrence_A(:), Coverage_A(:), ...
            Duration_B(:) * 1000, Occurrence_B(:), Coverage_B(:), ...
            Duration_C(:) * 1000, Occurrence_C(:), Coverage_C(:), ...
            Duration_D(:) * 1000, Occurrence_D(:), Coverage_D(:), ...
            Duration_E(:) * 1000, Occurrence_E(:), Coverage_E(:), ...
            Duration_F(:) * 1000, Occurrence_F(:), Coverage_F(:), ...
            Duration_G(:) * 1000, Occurrence_G(:), Coverage_G(:), ...
            'VariableNames', {'SubjectID', 'Group', 'Age', 'Gender', ...
            'Duration_A', 'Occurrence_A', 'Coverage_A', ...
            'Duration_B', 'Occurrence_B', 'Coverage_B', ...
            'Duration_C', 'Occurrence_C', 'Coverage_C', ...
            'Duration_D', 'Occurrence_D', 'Coverage_D', ...
            'Duration_E', 'Occurrence_E', 'Coverage_E', ...
            'Duration_F', 'Occurrence_F', 'Coverage_F', ...
            'Duration_G', 'Occurrence_G', 'Coverage_G'});

        T = rmmissing(T);

        % Ensure factors

        T.Group = removecats(categorical(T.Group));
        T.Gender = removecats(categorical(T.Gender)); % e.g., 'F' / 'M'

        % Grand-mean center Age (across all subjects)

        T.Age_c = T.Age - mean(T.Age, 'omitnan');

        if Sens_meds

            if strcmp(medSpec, 'classes')
                T.AD = AD';
                T.AP = AP';
                T.MS = MS';
                T.ANX = ANX';
                T.OTHER = OTHER';

                medNames = {'AD', 'AP', 'MS', 'ANX', 'OTHER'};
                medTerm = strjoin(medNames, ' + ');
            else
                T.poly_count = poly_count';
                medTerm = 'poly_count';
            end
        end

        metrics = {'Duration', 'Occurrence', 'Coverage'};
        letters = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};

        alphaFDR = 0.05; % FDR level for effects across the 7 microstates

        ResultsLME = struct();

        clc;

        for m = 1: numel(metrics)

            metric = metrics{m};

            fprintf('\n=============================\n');
            fprintf('Metric: %s (LMM per microstate)\n', metric);
            fprintf('=============================');

            lmePerState = cell(1, numel(letters));
            pGroup = nan(1, numel(letters));
            pAge = nan(1, numel(letters));
            pGender = nan(1, numel(letters));
            pAD = nan(1, numel(letters));
            pAP = nan(1, numel(letters));
            pMS = nan(1, numel(letters));
            pANX = nan(1, numel(letters));
            pOTHER = nan(1, numel(letters));
            pPOLY_COUNT = nan(1, numel(letters));

            % Fit 7 LMMs, one per microstate

            for s = 1: numel(letters)

                varName = sprintf('%s_%s', metric, letters{s});

                if ~Sens_meds
                    formula = sprintf('%s ~ Group + Age_c + Gender', varName);

                else
                    formula = sprintf('%s ~ Group + Age_c + Gender + %s', varName, medTerm);
                end

                % Fit LMM

                lme = fitlme(T, formula, 'FitMethod', 'REML');

                % Omnibus tests for Group, Age, Gender

                A = anova(lme, 'DFMethod', 'Satterthwaite');

                % Row names can vary a bit; find by name

                pGroup(s) = A.pValue(strcmp(string(A.Term), 'Group'));
                pAge(s) = A.pValue(strcmp(string(A.Term), 'Age_c'));
                pGender(s) = A.pValue(strcmp(string(A.Term), 'Gender'));

                if Sens_meds
                    if strcmp(medSpec, 'classes')
                        pAD(s) = A.pValue(strcmp(string(A.Term), 'AD'));
                        pAP(s) = A.pValue(strcmp(string(A.Term), 'AP'));
                        pMS(s) = A.pValue(strcmp(string(A.Term), 'MS'));
                        pANX(s) = A.pValue(strcmp(string(A.Term), 'ANX'));
                        pOTHER(s) = A.pValue(strcmp(string(A.Term), 'OTHER'));
                    else
                        pPOLY_COUNT(s) = A.pValue(strcmp(string(A.Term), 'poly_count'));
                    end

                end
                lmePerState{s} = lme;
            end

            % FDR across the 7 microstates for each effect

            alphaGroup = bh_fdr(pGroup, 0.05);
            alphaAge = bh_fdr(pAge, 0.05);
            alphaGender = bh_fdr(pGender, 0.05);

            if Sens_meds
                if strcmp(medSpec, 'classes')
                    alphaAD = bh_fdr(pAD, 0.05);
                    alphaAP = bh_fdr(pAP, 0.05);
                    alphaMS = bh_fdr(pMS, 0.05);
                    alphaANX = bh_fdr(pANX, 0.05);
                    alphaOTHER = bh_fdr(pOTHER, 0.05);
                else
                    alphaPOLY_COUNT = bh_fdr(pPOLY_COUNT, 0.05);
                end
            end

            % Report + posthocs (pairwise Group contrasts where alphaGroup < alpha)

            if ~Sens_meds
                report_LME_results_microstates(lmePerState, T, lower(metric), letters, alphaAge, alphaGender, alphaGroup, 0.05);
            else
                if strcmp(medSpec, 'classes')
                    report_LME_results_microstates_meds_classes(lmePerState, T, lower(metric), letters, alphaAge, alphaGender, alphaAD, alphaAP, alphaMS, alphaANX, alphaOTHER, alphaGroup, 0.05);
                else
                    report_LME_results_microstates_meds_poly(lmePerState, T, lower(metric), letters, alphaAge, alphaGender, alphaPOLY_COUNT, alphaGroup, 0.05);
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if test_transition_probability

        %% ======================================================================= %%
        %  TASK MICROSTATES — TWO TRANSITION ANALYSES (AGE-ADJUSTED, FDR CONTROL)
        %  Q1. Within each cell (Group), which X->Y transitions
        %         deviate from independence?  [Poisson GLMs with independence offset]
        %  Q2. Anchored contrasts: For each cell, which X->Y transitions in
        %         each target group differ from a chosen reference (age-anchored)?
        %
        %  Outputs:
        %   - TransStats_cellFDR: per-cell deviations (q-values, edges for plots)
        %   - CompStats: Target vs Reference anchored contrasts (q-values, edges)
        %%% ====================================================================== %%

        load('D:\Papers\2025\Submitted\Journal of Bipolar Disorders (resting state microstates)\Data\Participants_demographics.mat');
        load([saveDir '\outputStats.mat']); % per-subject MSClass etc.

        %% Test for an effect of group and microstate on mean microstate duration, mean microstate occurrence, and mean microstate coverage

        OccurrenceNames = {'Occurrence_A', 'Occurrence_B', 'Occurrence_C', 'Occurrence_D', 'Occurrence_E', 'Occurrence_F', 'Occurrence_G'};
        DurationNames = {'Duration_A', 'Duration_B', 'Duration_C', 'Duration_D', 'Duration_E', 'Duration_F', 'Duration_G'};

        % Collect age, gender and subject_ID vectors aligned to outputStats ordering

        AgeMap = struct(); GenderMap = struct; SubjectIDMap = struct;

        AgeMap.BP_I_Depressed = cellfun(@(s) s.age, France.BP_I.Depressed);
        GenderMap.BP_I_Depressed = cellfun(@(s) s.gender, France.BP_I.Depressed);
        SubjectIDMap.BP_I_Depressed = cellfun(@(s) s.number, France.BP_I.Depressed);

        AgeMap.BP_I_Euthymic = cellfun(@(s) s.age, France.BP_I.Euthymic);
        GenderMap.BP_I_Euthymic = cellfun(@(s) s.gender, France.BP_I.Euthymic);
        SubjectIDMap.BP_I_Euthymic = cellfun(@(s) s.number, France.BP_I.Euthymic);

        AgeMap.BP_II_Depressed = cellfun(@(s) s.age, France.BP_II.Depressed);
        GenderMap.BP_II_Depressed = cellfun(@(s) s.gender, France.BP_II.Depressed);
        SubjectIDMap.BP_II_Depressed = cellfun(@(s) s.number, France.BP_II.Depressed);

        AgeMap.BP_II_Euthymic = cellfun(@(s) s.age, France.BP_II.Euthymic);
        GenderMap.BP_II_Euthymic = cellfun(@(s) s.gender, France.BP_II.Euthymic);
        SubjectIDMap.BP_II_Euthymic = cellfun(@(s) s.number, France.BP_II.Euthymic);

        AgeMap.HC = cellfun(@(s) s.age, France.HC);
        GenderMap.HC = cellfun(@(s) s.gender, France.HC);
        SubjectIDMap.HC = cellfun(@(s) s.number, France.HC);

        AgeMap.Siblings = cellfun(@(s) s.age, France.Siblings);
        GenderMap.Siblings = cellfun(@(s) s.gender, France.Siblings);
        SubjectIDMap.Siblings = cellfun(@(s) s.number, France.Siblings);
  
        Out = struct();
    MicrostateLabels = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};

        for g = 1: nGroups
            
            grp = Groups{g};
            subj_list = outputStats.(grp);
            nSubj = numel(subj_list);

            for sj = 1: nSubj

                S = subj_list(sj);

                Transition_matrix = zeros(7, 7);
                Transition_matrix(1: 8: 49) = NaN;

                aa = [S.MSClass; -999]; % sentinel
                nonrepeats = find(diff(aa) ~= 0);
                segments = aa(nonrepeats).';
                segments = segments(segments ~= 0); % number of labels in a trial
                N_transitions = max(0, numel(segments) - 1);

                for p = 1:(numel(segments) - 1)
                    i = segments(p);
                    j = segments(p + 1);

                    if i ~= j
                        Transition_matrix(i, j) = Transition_matrix(i, j) + 1;
                    end
                end

                M_segments = numel(segments); % total number of labels

                C_segments = zeros(1, 7);

                for s = 1: 7
                    C_segments(s) = nnz(segments == s); % occurrence count of each label
                end

                Out.(grp){sj}.trans_counts = Transition_matrix;
                Out.(grp){sj}.M_segments = M_segments;
                Out.(grp){sj}.N_transitions = N_transitions;
                Out.(grp){sj}.C_segments = C_segments;
                Out.(grp){sj}.Duration_A = outputStats.(grp)(sj).MeanDuration(1);

                for s = 1: 7
                    Out.(grp){sj}.(DurationNames{s}) = outputStats.(grp)(sj).MeanDuration(s);
                    Out.(grp){sj}.(OccurrenceNames{s}) = outputStats.(grp)(sj).MeanOccurrence(s);
                end

                Out.(grp){sj}.transition_count_matrix = Transition_matrix / max(1, N_transitions);
            end
        end

        [srcIdx, tgtIdx] = find(~eye(7)); % 42 directed pairs
        GLM_OPTS = statset('MaxIter', 5000, 'TolX', 1e-10, 'TolFun', 1e-10);

        % Numerical guards (used below)

        LOGE_CLAMP = 12; % clamp logE to [-12, +12] to stabilize IRLS
        MIN_OBS = 3; % min pooled observed counts to fit a pair
        MIN_EXP = 3; % min pooled expected (or anchored mean) to fit

        % ---------------------------------------------------------------------------------------------------------------------------------------------------
        % Poisson GLMM with offset for microstate transitions, per group and per transition, adjusted for age and gender.
        % ---------------------------------------------------------------------------------------------------------------------------------------------------

        TransStats_cellFDR = struct();

        for g = 1: nGroups
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

            rows_Gender = [];

            for pop = 1: length(GenderMap.(grp))
                if strcmp(GenderMap.(grp)(pop), 'F')
                    rows_Gender = [rows_Gender; 0];
                else
                    rows_Gender = [rows_Gender; 1];
                end
            end

            cellOut = Out.(grp);

            IRR = NaN(7);
            beta0 = NaN(7);
            se = NaN(7);
            pval = NaN(7);
            ciLo = NaN(7);
            ciHi = NaN(7);
            nRows = zeros(7);
            ObsOverExp = NaN(7);

            RowsAll = table(); % y, logE, AgeZ, Gender, src, tgt

            if ~isempty(cellOut)
                for sj = 1: numel(cellOut)

                    S = cellOut{sj};
                    P = (S.C_segments(:) / S.M_segments);

                    Ssq = sum(P.^2);
                    den = 1 - Ssq;
                    kK = numel(S.C_segments);

                    gamma = 1 / den;

                    E = gamma * S.N_transitions * (P * P.');
                    E(1: kK + 1: kK  *kK) = NaN;

                    yvec = S.trans_counts(sub2ind([kK, kK], srcIdx, tgtIdx));
                    evec = E(sub2ind([kK, kK], srcIdx, tgtIdx));
                    keep = isfinite(evec) & evec > 0 & ~isnan(yvec);

                    if ~any(keep)
                        continue;
                    end

                    logE = log(evec(keep));
                    logE = max(min(logE, LOGE_CLAMP), -LOGE_CLAMP); % clamp

                    R = table(yvec(keep), logE, repmat(AgesZ(sj), nnz(keep), 1), repmat(rows_Gender(sj), nnz(keep), 1), srcIdx(keep), tgtIdx(keep), 'VariableNames', {'y', 'logE', 'AgeZ', 'Gender', 'src', 'tgt', 'Dur_src', 'Dur_trg'});
                    RowsAll = [RowsAll; R];
                end
            end

            for k = 1: numel(srcIdx)
                s = srcIdx(k);
                t = tgtIdx(k);
                R = RowsAll(RowsAll.src == s & RowsAll.tgt == t, :);

                % screens

                if sum(R.y) < MIN_OBS || sum(exp(R.logE)) < MIN_EXP
                    continue;
                end

                nRows(s, t) = height(R);
                ObsOverExp(s, t) = sum(R.y) / sum(exp(R.logE));

                try
                    mdl = fitglm(R, 'y ~ 1 + AgeZ + Gender', 'Distribution', 'poisson', 'Link', 'log', 'Offset', R.logE, 'Options', GLM_OPTS);

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

                catch
                end
            end

            [sigMask, qMat] = bh_fdr_mask(pval, 0.05);

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
                edges.(grp) = table(sSig, tSig, irr_vec, lo_vec, hi_vec, o2e_vec, q_vec, 'VariableNames',{'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'ObsOverExp', 'qFDR'});
                edges.(grp).srcLabel = reshape(string(MicrostateLabels(edges.(grp).src)), [], 1);
                edges.(grp).tgtLabel = reshape(string(MicrostateLabels(edges.(grp).tgt)), [], 1);
                dlab = strings(numel(sSig), 1);
                dlab(irr_vec > 1) = "Above expected";
                dlab(irr_vec < 1) = "Below expected";
                edges.(grp).Direction = dlab;
            else
                edges.(grp) = table([], [], [], [], [], [], [], 'VariableNames',{'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'ObsOverExp', 'qFDR'});
                edges.(grp).srcLabel = strings(0, 1);
                edges.(grp).tgtLabel = strings(0, 1);
                edges.(grp).Direction = strings(0, 1);
            end

            TransStats_cellFDR.(grp).source = srcIdx;
            TransStats_cellFDR.(grp).target = tgtIdx;
            TransStats_cellFDR.(grp).IRR = IRR;
            TransStats_cellFDR.(grp).beta0 = beta0;
            TransStats_cellFDR.(grp).se = se;
            TransStats_cellFDR.(grp).p = pval;
            TransStats_cellFDR.(grp).q = qMat;
            TransStats_cellFDR.(grp).sig = sigMask;
            TransStats_cellFDR.(grp).dir = dirMatrix;
            TransStats_cellFDR.(grp).nRows = nRows;
            TransStats_cellFDR.(grp).ObsOverExp = ObsOverExp;
            TransStats_cellFDR.(grp).edges = edges.(grp);

        end

        clear Tstats Tstats_pres

        % Plot p-values as heatmaps

        min_pFDR = min([TransStats_cellFDR.BP_I_Depressed.q(:); TransStats_cellFDR.BP_I_Euthymic.q(:); TransStats_cellFDR.BP_II_Depressed.q(:); TransStats_cellFDR.BP_II_Euthymic.q(:); TransStats_cellFDR.HC.q(:); TransStats_cellFDR.Siblings.q(:)]);
        max_pFDR = max([TransStats_cellFDR.BP_I_Depressed.q(:); TransStats_cellFDR.BP_I_Euthymic.q(:); TransStats_cellFDR.BP_II_Depressed.q(:); TransStats_cellFDR.BP_II_Euthymic.q(:); TransStats_cellFDR.HC.q(:); TransStats_cellFDR.Siblings.q(:)]);

        figure('Name', 'BH-FDR Heatmaps: all groups', 'NumberTitle', 'off', 'Color', [1 1 1]);
        tiled = tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

        for g = 1: numel(Groups)

            grp = Groups{g};

            % Create an RGB image for the heatmap

            RGB = NaN(7, 7, 3);
            whiteColor = [1 1 1];

            for source = 1: 7
                for target = 1: 7
                    if source == target
                        RGB(source, target, :) = whiteColor; % diagonal
                    else
                        pval = TransStats_cellFDR.(grp).q(source, target);

                        if isnan(pval) || pval > 0.05

                            % Non-significant: blue color (gradient by p)

                            RGB(source, target, :) = blueColor(pval, min_pFDR, max_pFDR); % define below
                        else

                            % Significant: red color (gradient by p)

                            RGB(source, target, :) = redColor(pval, min_pFDR, max_pFDR); % define below
                        end
                    end
                end
            end

            % Plot using image (exact grid)

            ax = nexttile(g, [1 1]);
            imshow(RGB, 'InitialMagnification', 'fit', 'Parent', ax);
            axis(ax, 'on'); axis(ax, 'image');

            % Axes labels

            set(ax, 'XTick', 1: 7, 'XTickLabel', MicrostateLabels, 'YTick', 1: 7, 'YTickLabel', MicrostateLabels, 'Box', 'on', 'FontName', 'TimesNewRoman', 'FontSize', 13, 'FontWeight', 'normal');
            ylabel('Origin', 'FontName', 'TimesNewRoman', 'FontSize', 14, 'FontWeight', 'bold'); xlabel('Target', 'FontName', 'TimesNewRoman', 'FontSize', 14, 'FontWeight', 'bold');

            % Draw thin black grid lines and thick frames for significant cells

            hold(ax, 'on');

            % Grid lines (subtle)

            for k = 0.5: 1: 6.5
                plot(ax, [k k], [0.5 7.5], 'k-', 'LineWidth', 1);
                plot(ax, [0.5 7.5], [k k], 'k-', 'LineWidth', 1);
            end

            % Thick black frames around significant cells

            for source = 1: 7
                for target = 1: 7

                    if i == j
                        continue;
                    end

                    if TransStats_cellFDR.(grp).q(source, target) <= 0.05 && ~isnan(TransStats_cellFDR.(grp).q(source, target))

                        % Draw a rectangle around cell (i,j)

                        x = target - 0.5; y = source - 0.5;
                        w = 1; h = 1;
                        rectangle(ax, 'Position', [x, y, w, h], 'EdgeColor', 'k', 'LineWidth', 4);
                    end
                end
            end

            hold(ax, 'off');
        end

        Ngrad = 256;

        pBlue = linspace(0.05, max_pFDR, Ngrad/2);
        gradBlue = zeros(length(pBlue), 1, 3);

        for k = 1: length(pBlue)
            gradBlue(k, 1, :) = blueColor(pBlue(k), min_pFDR, max_pFDR);
        end

        % Build red ramp (bottom: min_pFDR, top: p = 0.05)

        pRed = linspace(min_pFDR, 0.05, Ngrad/2);
        gradRed = zeros(length(pRed), 1, 3);

        for k = 1:length(pRed)
            gradRed(k, 1, :) = redColor(pRed(k), min_pFDR, max_pFDR);
        end

        % Combine into a single gradient: blue (bottom) + red (top)

        gradFull = cat(1, gradRed, gradBlue);
        gradFull = max(0, min(1, gradFull)); % ensure in [0,1]

        % Create a single colorbar as an image (no gaps)

        cbAx = axes('Position', [0.96, 0.802, 0.02, 0.16]);
        cbImg = reshape(gradFull, [size(gradFull,1), 1, 3]);
        imagesc(cbImg, 'Parent', cbAx);
        set(cbAx, 'YDir', 'normal');
        axis(cbAx, 'off');

        text(1.8, 0.5, sprintf('%.3f', min_pFDR), 'FontName', 'TimesNewRoman', 'FontSize', 12); text(1.8, 128, num2str(0.05), 'FontName', 'TimesNewRoman', 'FontSize', 12); text(1.8, 256, sprintf('%.3f', max_pFDR), 'FontName', 'TimesNewRoman', 'FontSize', 12);

        %%%% Plot a circular plot that shows that transitions' p value, polarity (above or below expeced), and duration spent in sounrce and target microstates

        vals_dur = [];

        try

            clear vals1 vals2

            vals1 = cellfun(@(s) s.duration_start(:), transitions_Siblings, 'UniformOutput', false);
            vals2 = cellfun(@(s) s.duration_end(:), transitions_Siblings, 'UniformOutput', false);

            for kkk = 1: size(vals1, 2)
                vals_dur = [vals_dur; nanmean(vals1{kkk}); nanmean(vals2{kkk})];
            end

        catch
        end


    end

end

function Check_demographics(France)

% Look for normality of age and scores using Shpairo-Wilk test, used for small sample sizes (< 50 samples) and
% check for significant difference across groups using ANOVA1 or its non parametric equivalent - kruskalwallis.

index = 0;

for n = 1: length(France.BP_I.Depressed)
    index = index + 1;
    Age(index) = France.BP_I.Depressed{n}.age;
    Gender(index) = France.BP_I.Depressed{n}.gender;
    GAF(index) = France.BP_I.Depressed{n}.GAF;
    MADRS(index) = France.BP_I.Depressed{n}.MADRS;
    YMRS(index) = France.BP_I.Depressed{n}.YMRS;

    Group{index} = 'BP_I_Depressed';

end

index1 = index;

for n = 1: length(France.BP_I.Euthymic)
    index1 = index1 + 1;
    Age(index1) = France.BP_I.Euthymic{n}.age;
    Gender(index1) = France.BP_I.Euthymic{n}.gender;
    GAF(index1) = France.BP_I.Euthymic{n}.GAF;
    MADRS(index1) = France.BP_I.Euthymic{n}.MADRS;
    YMRS(index1) = France.BP_I.Euthymic{n}.YMRS;

    Group{index1} = 'BP_I_Euthymic';

end

index2 = index1;

for n = 1: length(France.BP_II.Depressed)
    index2 = index2 + 1;
    Age(index2) = France.BP_II.Depressed{n}.age;
    Gender(index2) = France.BP_II.Depressed{n}.gender;
    GAF(index2) = France.BP_II.Depressed{n}.GAF;
    MADRS(index2) = France.BP_II.Depressed{n}.MADRS;
    YMRS(index2) = France.BP_II.Depressed{n}.YMRS;

    Group{index2} = 'BP_II_Depressed';

end

index3 = index2;

for n =1: length(France.BP_II.Euthymic)
    index3 = index3 + 1;
    Age(index3) = France.BP_II.Euthymic{n}.age;
    Gender(index3) = France.BP_II.Euthymic{n}.gender;
    GAF(index3) = France.BP_II.Euthymic{n}.GAF;
    MADRS(index3) = France.BP_II.Euthymic{n}.MADRS;
    YMRS(index3) = France.BP_II.Euthymic{n}.YMRS;

    Group{index3} = 'BP_II_Euthymic';

end

index4 = index3;

for n = 1: length(France.HC)
    index4 = index4 + 1;
    Age(index4) = France.HC{n}.age;
    Gender(index4) = France.HC{n}.gender;
    GAF(index4) = 0/0;
    MADRS(index4) = 0/0;
    YMRS(index4) = 0/0;

    Group{index4} = 'HC';

end

index5 = index4;

for n = 1: length(France.Siblings)
    index5 = index5 + 1;
    Age(index5) = France.Siblings{n}.age;
    Gender(index5) = France.Siblings{n}.gender;
    GAF(index5) = 0/0;
    MADRS(index5) = 0/0;
    YMRS(index5) = 0/0;

    Group{index5} = 'Siblings';

end

disp('------------------------------------------- Demographics ------------------------------------------')
disp(' ');

disp(' _________________Age___________________');
disp(' ');
disp(['BP_I Depressed: ' num2str(round(mean(Age(1: index)) * 10) / 10) ' +- ' num2str(round(std(Age(1: index)) * 10) / 10)]);
disp(' ');
disp('Is age normally distributed?')
disp(' ');
normalitytest(Age(1: index));
disp(' ');
disp(' ');
disp(' ');

disp(['BP_I Euthymic: ' num2str(round(mean(Age(index + 1: index1)) * 10) / 10) ' +- ' num2str(round(std(Age(index + 1: index1)) * 10) / 10)]);
disp(' ');
disp('Is age normally distributed?')
disp(' ');
normalitytest(Age(index + 1: index1));
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Depressed: ' num2str(round(mean(Age(index1 + 1: index2)) * 10) / 10) ' +- ' num2str(round(std(Age(index1 + 1: index2)) * 10) / 10)]);
disp(' ');
disp('Is age normally distributed?')
disp(' ');
normalitytest(Age(index1 + 1: index2));
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Euthymic: '  num2str(round(mean(Age(index2 + 1: index3)) * 10) / 10) ' +- ' num2str(round(std(Age(index2 + 1: index3)) * 10) / 10)]);
disp(' ');
disp('Is age normally distributed?')
disp(' ');
normalitytest(Age(index2 + 1: index3));
disp(' ');
disp(' ');
disp(' ');

disp(['HC: ' num2str(round(mean(Age(index3 + 1: index4)) * 10) / 10) ' +- ' num2str(round(std(Age(index3 + 1: index4)) * 10) / 10)]);
disp(' ');
disp('Is age normally distributed?')
disp(' ');
normalitytest(Age(index3 + 1: index4));
disp(' ');
disp(' ');
disp(' ');

disp(['Siblings: ' num2str(round(mean(Age(index4 + 1: index5)) * 10) / 10) ' +- ' num2str(round(std(Age(index4 + 1: index5)) * 10) / 10)]);
disp(' ');
disp('Is age normally distributed?')
disp(' ');
Results = normalitytest(Age(index4 + 1: index5));
disp(' ');
disp(' ');
disp(' ');

[~, ~, stats] = anova1(Age, Group);
[results, ~, ~, gnames] = multcompare(stats);
tbl = array2table(results, "VariableNames", ["Group", "Control Group", "Lower Limit", "Difference", "Upper Limit", "P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

disp(' _________________Gender ___________________');
disp(' ');
disp(['BP_I Depressed: M - ' num2str(length(find(Gender(1: index) == 'M'))) ', F - '  num2str(length(find(Gender(1: index) == 'F')))]);
disp(['BP_I Euthymic: M - ' num2str(length(find(Gender(index + 1: index1) == 'M'))) ', F - '  num2str(length(find(Gender(index + 1: index1) == 'F')))]);
disp(['BP_II Depressed: M - ' num2str(length(find(Gender(index1 + 1: index2) == 'M'))) ', F - '  num2str(length(find(Gender(index1 + 1: index2) == 'F')))]);
disp(['BP_II Euthymic: M - ' num2str(length(find(Gender(index2 + 1: index3) == 'M'))) ', F - '  num2str(length(find(Gender(index2 + 1: index3) == 'F')))]);
disp(['HC: M - ' num2str(length(find(Gender(index3 + 1: index4) == 'M'))) ', F - '  num2str(length(find(Gender(index3 + 1: index4) == 'F')))]);
disp(['Siblings: M - ' num2str(length(find(Gender(index4 + 1: index5) == 'M'))) ', F - '  num2str(length(find(Gender(index4 + 1: index5) == 'F')))]);

disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(' _________________GAF ___________________');
disp(' ');
disp(['BP_I Depressed: ' num2str(round(mean(GAF(1: index)) * 10) / 10) ' +- ' num2str(round(std(GAF(1: index)) * 10) / 10)]);
disp(' ');
disp('Is GAF normally distributed?')
disp(' ');
clear aaa;
aaa = GAF(1: index);
aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_I Euthymic: ' num2str(round(mean(GAF(index + 1: index1)) * 10) / 10) ' +- ' num2str(round(std(GAF(index + 1: index1)) * 10) / 10)]);
disp(' ');
disp('Is GAF normally distributed?')
disp(' ');
clear aaa;
aaa = GAF(index + 1: index1);
aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Depressed: ' num2str(round(mean(GAF(index1 + 1: index2)) * 10) / 10) ' +- ' num2str(round(std(GAF(index1 + 1: index2)) * 10) / 10)]);
disp(' ');
disp('Is GAF normally distributed?')
disp(' ');
clear aaa;
aaa = GAF(index1 + 1: index2);
aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Euthymic: '  num2str(round(mean(GAF(index2 + 1: index3)) * 10) / 10) ' +- ' num2str(round(std(GAF(index2 + 1: index3)) * 10) / 10)]);
disp(' ');
disp('Is GAF normally distributed?')
disp(' ');
clear aaa;
aaa = GAF(index2 + 1: index3);
aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

clear stats results gnames

[~, ~, stats] = anova1(GAF, Group);
[results, ~, ~, gnames] = multcompare(stats);
tbl = array2table(results,"VariableNames", ["Group", "Control Group", "Lower Limit", "Difference", "Upper Limit", "P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(' _________________MADRS ___________________');
disp(' ');
disp(['BP_I Depressed: ' num2str(round(mean(MADRS(1: index)) * 10) / 10) ' +- ' num2str(round(std(MADRS(1: index)) * 10) / 10)]);
disp(' ');
disp('Is MADRS normally distributed?')
disp(' ');
clear aaa;
aaa = MADRS(1: index);
aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_I Euthymic: ' num2str(round(mean(MADRS(index + 1: index1)) * 10) / 10) ' +- ' num2str(round(std(MADRS(index + 1: index1)) * 10) / 10)]);
disp(' ');
disp('Is MADRS normally distributed?')
disp(' ');
clear aaa;
aaa = MADRS(index + 1: index1);
aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Depressed: ' num2str(round(mean(MADRS(index1 + 1: index2)) * 10) / 10) ' +- ' num2str(round(std(MADRS(index1 + 1: index2)) * 10) / 10)]);
disp('Is MADRS normally distributed?')
disp(' ');
clear aaa;
aaa = MADRS(index1 + 1: index2);
aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Euthymic: '  num2str(round(mean(MADRS(index2 + 1: index3)) * 10) / 10) ' +- ' num2str(round(std(MADRS(index2 + 1: index3)) * 10) / 10)]);
disp('Is MADRS normally distributed?')
disp(' ');
clear aaa;
aaa = MADRS(index2 + 1: index3);
aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

clear stats results gnames

[~, ~, stats] = anova1(MADRS, Group);
[results, ~, ~, gnames] = multcompare(stats);
tbl = array2table(results,"VariableNames", ["Group", "Control Group", "Lower Limit", "Difference", "Upper Limit", "P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(' _________________YMRS ___________________');
disp(' ');
disp(['BP_I Depressed: ' num2str(round(mean(YMRS(1: index)) * 10) / 10) ' +- ' num2str(round(std(YMRS(1: index)) * 10) / 10)]);
disp('Is YMRS normally distributed?')
disp(' ');
clear aaa;
aaa = YMRS(1: index);
aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_I Euthymic: ' num2str(round(mean(YMRS(index + 1: index1)) * 10) / 10) ' +- ' num2str(round(std(YMRS(index + 1: index1)) * 10) / 10)]);
disp('Is YMRS normally distributed?')
disp(' ');
clear aaa;
aaa = YMRS(index + 1: index1);
aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Depressed: ' num2str(round(mean(YMRS(index1 + 1: index2)) * 10) / 10) ' +- ' num2str(round(std(YMRS(index1 + 1: index2)) * 10) / 10)]);
disp('Is YMRS normally distributed?')
disp(' ');
clear aaa;
aaa = YMRS(index1 + 1: index2);
aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Euthymic: '  num2str(round(mean(YMRS(index2 + 1: index3)) * 10) / 10) ' +- ' num2str(round(std(YMRS(index2 + 1: index3)) * 10) / 10)]);
disp('Is YMRS normally distributed?')
disp(' ');
clear aaa;
aaa = YMRS(index2 + 1: index3);
aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

clear stats results gnames

[~, ~, stats] = anova1(YMRS, Group);
[results, ~, ~, gnames] = multcompare(stats);
tbl = array2table(results, "VariableNames", ["Group", "Control Group", "Lower Limit", "Difference", "Upper Limit", "P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

end

function c = blueColor(p, min_pFDR, max_pFDR)

p_min_blue = 0.05;
p_max_blue = max_pFDR;

p_clamped = min(max(p, p_min_blue), p_max_blue);

if p_max_blue > p_min_blue
    t = (p_clamped - p_min_blue) / (p_max_blue - p_min_blue);
else
    t = 0;
end

paleBlue = [0.95, 0.98, 1.00];
darkBlue = [0.00, 0.20, 0.70];

c = (1 - t) * paleBlue + t * darkBlue;
c = max(0, min(1, c));

end

function c = redColor(p, min_pFDR, max_pFDR)

p_min_signif = min_pFDR;
p_max_signif = 0.05;

p_clamped = min(max(p, p_min_signif), p_max_signif);

if p_max_signif > p_min_signif
    t = (p_max_signif - p_clamped) / (p_max_signif - p_min_signif);
else
    t = 0;
end

paleRed = [1.0, 0.92, 0.92];
darkRed = [1.0, 0.00, 0.00];

c = (1 - t) * paleRed + t * darkRed;
c = max(0, min(1, c));
end

function Compute_and_plot_individ_total_GEV(GEV)

if ~exist('GEV', 'var')
    load('D:\Papers\2025\Submitted\Journal of Bipolar Disorders (resting state microstates)\Analysis\Individual_micro_states_3_15.mat');
end

% Your preferred display order (only plot those that exist)

preferredOrder = {'BP_I_Depressed', 'BP_I_Euthymic', 'BP_II_Depressed', 'BP_II_Euthymic', 'HC', 'Siblings'};
showOrder = {};

for g = 1: numel(preferredOrder)
    if isfield(GEV, preferredOrder{g})
        showOrder{end + 1} = preferredOrder{g};
    end
end

figure(1); set(gcf, 'Color', [1 1 1]);

ymax = 0.75; ymin = 0.35;  % keep your axes

for i = 1: min(6, numel(showOrder))
    fld = showOrder{i};
    mu  = MeanGEV.(fld);
    sd  = StdGEV.(fld);

    % Print percent changes (K-1 vs K)

    for n = (Ks(1) + 1): Ks(end)
        pct = ((mu(n - 2) - mu(n - 3)) / mu(n - 3)) * 100;
        fprintf('%s. %d vs. %d: %.1f%%\n', strrep(fld, '_', '\_'), n - 1, n, pct);
    end

    fprintf('\n\n');

    subplot(3, 2, i);

    try
        shadedErrorBar(3: length(mu) + 2, mu, sd);
    catch
        plot(3: length(mu) + 2, mu, 'LineWidth', 1.5); hold on;
        plot(3: length(mu) + 2, mu + sd, '--'); plot(3: length(mu) + 2, mu - sd, '--'); hold off;
    end
    xlim([Ks(1) Ks(end)]);
    ylim([ymin ymax]); title(strrep(fld, '_', '\_'));
end

% "7 microstate templates" summary (if available)

if any(Ks == 7)

    k7 = 7;

    for i = 1:numel(showOrder)
        fld = showOrder{i};
        vals = GEV.(fld)(:, k7);
        fprintf('7 microstate templates: %s GEV = %d%% %s %d%%\n', strrep(fld, '_', '\_'), round(nanmean(vals) * 100), char(177), round(nanstd(vals) * 100));
    end
end
end

function q = bh_fdr(p, alpha)

p = p(:);
n = numel(p);
[ps, idx] = sort(p);
ranks = (1: n)';
qtemp = (ps .* n) ./ ranks;
qtemp = cummin(flipud(qtemp));
qtemp = flipud(qtemp);
q = nan(size(p));
q(idx) = min(qtemp, 1);

if nargin > 1
end

end
