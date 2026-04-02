function ERP_stats = ERP_validation_stats(dataPath)

% ERP validation stats: Subject-level ROI mean amplitudes in predefined windows.
% Groups are pooled to match Fig. 2: HC, BD, Siblings.

S = load(dataPath);
France = S.France;
clear S

conds = {'Negative', 'Neutral', 'Positive'};
groupNames = {'HC', 'BD', 'Siblings'};
groupCells = {France.HC, [France.BP_I.Depressed, France.BP_I.Euthymic, France.BP_II.Depressed, France.BP_II.Euthymic], France.Siblings};

winNames = {'N200', 'P300', 'LPP'};
winEdges = {[180 300], [300 500], [500 1000]};

roi.N200 = {'Fz', 'FC1', 'FC2', 'FC3', 'FC4', 'Cz', 'AFF1h', 'AFF2h'};
roi.P300 = {'CPz', 'CP1', 'CP2', 'Pz', 'P1', 'P2', 'POz', 'PO3', 'PO4'};
roi.LPP = roi.P300;

EEG0 = France.HC{1}.Neutral.EEG_Data;
chanlocs = EEG0.chanlocs;
times = EEG0.times;

baselineIdx = times >= -250 & times < 0;

roiIdx.N200 = local_chan_idx(chanlocs, roi.N200);
roiIdx.P300 = local_chan_idx(chanlocs, roi.P300);
roiIdx.LPP = local_chan_idx(chanlocs, roi.LPP);

Subject = {};
Group = {};
Condition = {};
Window = {};
Age = [];
Amplitude = [];

subCounter = 0;

for g = 1: numel(groupNames)

    subjCell = groupCells{g};

    for s = 1: numel(subjCell)

        subCounter = subCounter + 1;
        subj = subjCell{s};

        for c = 1: numel(conds)

            condName = conds{c};
            EEG = subj.(condName).EEG_Data;
            data = double(EEG.data);

            base = mean(data(:, baselineIdx, :), 2);
            data = data - base;

            for w = 1: numel(winNames)

                winName = winNames{w};
                edges = winEdges{w};
                winIdx = times >= edges(1) & times < edges(2);

                amp = mean_window_roi(data, roiIdx.(winName), winIdx);

                Subject{end + 1, 1} = sprintf('S%03d', subCounter);
                Group{end + 1, 1} = groupNames{g};
                Condition{end + 1, 1} = condName;
                Window{end + 1, 1} = winName;
                Age(end + 1, 1) = subj.age;
                Amplitude(end + 1, 1) = amp;
            end
        end
    end
end

T = table(categorical(Subject), categorical(Group), categorical(Condition), categorical(Window), Age, Amplitude, 'VariableNames', {'Subject', 'Group', 'Condition', 'Window', 'Age', 'Amplitude'});

ERP_stats = struct();

for w = 1: numel(winNames)

    thisWin = winNames{w};
    Tsub = T(T.Window == categorical(string(thisWin)), :);
    Tsub.Age_c = Tsub.Age - mean(Tsub.Age, 'omitnan');

    lme = fitlme(Tsub, 'Amplitude ~ Group * Condition + Age_c + (1 | Subject)');
    a = anova(lme, 'DFMethod', 'Satterthwaite');

    ERP_stats.(thisWin).table = Tsub;
    ERP_stats.(thisWin).lme = lme;
    ERP_stats.(thisWin).anova = a;

    fprintf('\n================ %s ================\n', thisWin);
    disp(a)
end

end

function idx = local_chan_idx(chanlocs, labels)

allLabs = string({chanlocs.labels});
idx = zeros(1, numel(labels));

for i = 1: numel(labels)

    hit = find(strcmpi(allLabs, labels{i}), 1, 'first');

    if isempty(hit)

        error('Channel label not found: %s', labels{i});

    end

    idx(i) = hit;
end

end

function amp = mean_window_roi(data, roiIdx, winIdx)

roiData = data(roiIdx, winIdx, :);
roiData = squeeze(mean(roiData, 1));

if isvector(roiData)

    amp = mean(roiData, 'omitnan');

else

    roiByTrial = mean(roiData, 1, 'omitnan');
    amp = mean(roiByTrial, 'omitnan');

end

end