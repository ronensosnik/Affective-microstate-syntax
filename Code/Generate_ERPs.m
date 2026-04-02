function Generate_ERPs(cfg)

cfg = set_defaults(cfg);

S = load(cfg.dataPath);
France = S.France;
clear S

S = load(cfg.outputStatsPath, 'outputStats');
outputStats = S.outputStats;
clear S

[chanlocs, timesAll] = get_chanlocs_times(France, cfg.conds{2});

roiIdxN200 = local_chan_idx(chanlocs, cfg.roiN200);
roiIdxP300 = local_chan_idx(chanlocs, cfg.roiP300_LPP);

baselineIdx = timesAll >= cfg.baselineMs(1) & timesAll < cfg.baselineMs(2);
plotIdx = timesAll >= cfg.plotMs(1) & timesAll <= cfg.plotMs(2);
tPlot = timesAll(plotIdx);

winIdxN200 = timesAll >= cfg.winN200(1) & timesAll < cfg.winN200(2);
winIdxP300 = timesAll >= cfg.winP300(1) & timesAll < cfg.winP300(2);
winIdxLPP = timesAll >= cfg.winLPP(1)  & timesAll < cfg.winLPP(2);

% ---- Build pooled-group subject lists for ERP/topographies ----

subjHC = France.HC;
subjSib = France.Siblings;
subjBD = [France.BP_I.Depressed, France.BP_I.Euthymic, France.BP_II.Depressed, France.BP_II.Euthymic];

subjAll = [subjHC, subjSib, subjBD];

groupNames = {'HC', 'BD', 'Siblings'};
groupCells = {subjHC, subjBD, subjSib};

% ---- ERP waveforms (per condition, per ROI, per group) ----

erpN200 = struct();
erpP300 = struct();

for c = 1: numel(cfg.conds)
    condName = cfg.conds{c};
    erpN200.(condName) = compute_erp_by_group(groupCells, condName, roiIdxN200, baselineIdx, plotIdx);
    erpP300.(condName) = compute_erp_by_group(groupCells, condName, roiIdxP300, baselineIdx, plotIdx);
end

% ==========================================================================================
% Peak topographies per ERP (N200 / P300), group (HC/BD/Siblings) and condition (Negative, Positive, Neutral).
% For each condition and group: find N200 peak latency (minimum of the group-mean FC ROI waveform) and P300 peak
% latency (maximum of the group-mean ROI), then plot the scalp voltage distribution at that time point.
% ==========================================================================================

% =======  N200 ==========

topoN200_peak = struct(); % topoN200_peak.(condName){g} = [nChan x 1]
tN200_peak_ms = struct(); % tN200_peak_ms.(condName)(g) = scalar ms

idxN200_plot = tPlot >= cfg.winN200(1) & tPlot < cfg.winN200(2);
idxN200_plot_list = find(idxN200_plot);

for c = 1: numel(cfg.conds)

    condName = cfg.conds{c};
    topoN200_peak.(condName) = cell(1, 3);
    tN200_peak_ms.(condName) = nan(1, 3);

    for g = 1: 3

        X = erpN200.(condName).sub{g}; % [nSub x nTime]
        mu = mean(X, 1, 'omitnan'); % [1 x nTime]
        [~, iRel] = min(mu(idxN200_plot)); % peak within 180–300 ms
        iPlot = idxN200_plot_list(iRel); % index into tPlot
        tPeak = tPlot(iPlot); % ms
        tN200_peak_ms.(condName)(g) = tPeak;

        [~, iTimeAll] = min(abs(timesAll - tPeak)); % nearest sample in full time vector
        topoN200_peak.(condName){g} = compute_topo_at_sample(groupCells{g}, condName, baselineIdx, iTimeAll);
    end
end

climN200_peak = get_clim_peak_struct(topoN200_peak, cfg.conds);

% ======= P300 =========

topoP300_peak = struct(); % topoP300_peak.(condName){g} = [nChan x 1]
tP300_peak_ms = struct(); % tP300_peak_ms.(condName)(g) = scalar ms

idxP300_plot = tPlot >= cfg.winP300(1) & tPlot < cfg.winP300(2);
idxP300_plot_list = find(idxP300_plot);

for c = 1: numel(cfg.conds)

    condName = cfg.conds{c};
    topoP300_peak.(condName) = cell(1, 3);
    tP300_peak_ms.(condName) = nan(1, 3);

    for g = 1: 3

        X = erpP300.(condName).sub{g}; % [nSub x nTime]
        mu = mean(X, 1, 'omitnan'); % [1 x nTime]
        [~, iRel] = max(mu(idxP300_plot)); % peak within 300–500 ms
        iPlot = idxP300_plot_list(iRel); % index into tPlot
        tPeak = tPlot(iPlot); % ms
        tP300_peak_ms.(condName)(g) = tPeak;

        [~, iTimeAll] = min(abs(timesAll - tPeak)); % nearest sample in full time vector
        topoP300_peak.(condName){g} = compute_topo_at_sample(groupCells{g}, condName, baselineIdx, iTimeAll);
    end
end

climP300_peak = get_clim_peak_struct(topoP300_peak, cfg.conds);

% ====================================================================================
% LPP mean-amplitude topographies (500–1000 ms), one map per condition,
% averaged across all participants in all groups (subjAll).
% ====================================================================================

topoLPP_mean = struct();

for c = 1: numel(cfg.conds)
    condName = cfg.conds{c};
    topoLPP_mean.(condName) = compute_topo_mean_window(subjAll, condName, baselineIdx, winIdxLPP);
end

valsLPP = [];

for c = 1: numel(cfg.conds)
    condName = cfg.conds{c};
    valsLPP = [valsLPP; topoLPP_mean.(condName)(:)];
end

mxLPP = max(abs(valsLPP));
climLPP_mean = [-mxLPP mxLPP];

% Use one common scale across N200, P300, and LPP topoplots

clim = [min([climN200_peak climP300_peak climLPP_mean]) max([climN200_peak climP300_peak climLPP_mean])];

% ---- Microstate label probability (Panel B) ----

msProb = compute_ms_prob(outputStats, cfg.msGroups, cfg.condsLower, cfg.plotMs);

% ---- Heatmap color limits (shared across all 9 panels) ----

allP = [];

for r = 1: numel(cfg.msGroups)
    grpName = cfg.msGroups{r};

    for c = 1: numel(cfg.condsLower)
        condName = cfg.condsLower{c};
        allP = [allP; msProb.(grpName).(condName)(:)];
    end
end

pHi = prctile(allP, 95);
pHi = max(pHi, 0.05);
heatCLim = [0 pHi];

%% -------------------------------------- Plot figure ---------------------------------------------

figure1 = figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05 0.05 0.90 0.90], 'InvertHardcopy', 'off');

x = 0.04;
w = 0.92;

yTop = 0.96;
yBottom = 0.03;
gapAB = 0.02;

totalH = yTop - yBottom - gapAB;

hA = totalH / 2;
hB = totalH / 2;

yA = yBottom + hB + gapAB;
yB = yBottom;

pA = uipanel(figure1, 'Units', 'normalized', 'Position', [x yA w hA], 'BorderType', 'none', 'BackgroundColor', 'w');
pB = uipanel(figure1, 'Units', 'normalized', 'Position', [x yB w hB], 'BorderType', 'none', 'BackgroundColor', 'w');

% ============== Panel A layout ==============

tA = tiledlayout(pA, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

figure(2);
t2 = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

figure(3);
t3 = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% LPP topographies in a separate figure, one per condition (all subjects)

figure(4);
t4 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for c = 1: 3

    condName = cfg.conds{c};

    ax1 = nexttile(tA, c);
    plot_erp(ax1, tPlot, erpN200.(condName), cfg.groupColors);
    grid('off');
    title(ax1, condName);
    ylabel(ax1, 'FC ROI (\muV)');

    ax2 = nexttile(tA, 3 + c);
    plot_erp(ax2, tPlot, erpP300.(condName), cfg.groupColors);
    grid('off');
    ylabel(ax2, 'CP–P ROI (\muV)');

    % N200 topography

    for g = 1: 3
        ttl = sprintf('%s %s (%d ms)', groupNames{g}, condName, round(tN200_peak_ms.(condName)(g)));
        args = {'electrodes', 'on', 'emarker', {'.', 'k', 10, 1}, 'emarker2', {roiIdxN200, '.', 'k', 20, 2}};
        ax3 = nexttile(t2, (g * 3) - 2 + (c - 1));
        axes(ax3);
        topoplot(topoN200_peak.(condName){g}, chanlocs, 'maplimits', clim, args{:});
        title(ax3, ttl)
    end

    % P300 topographies

    for g = 1: 3
        ttl = sprintf('%s %s (%d ms)', groupNames{g}, condName, round(tP300_peak_ms.(condName)(g)));
        args = {'electrodes', 'on', 'emarker', {'.', 'k', 10, 1}, 'emarker2', {roiIdxP300, '.', 'k', 20, 2}};
        ax4 = nexttile(t3, (g * 3) - 2 + (c - 1));
        axes(ax4);
        topoplot(topoP300_peak.(condName){g}, chanlocs, 'maplimits', clim, args{:});
        title(ax4, ttl)
    end

    % LPP mean-amplitude topography (all participants, all groups)

    ttl = sprintf('All groups %s (LPP 500 - 1000 ms)', condName);
    args = {'electrodes', 'on', 'emarker', {'.', 'k', 10, 1}, 'emarker2', {roiIdxP300, '.', 'k', 20, 2}};
    ax5 = nexttile(t4, c);
    axes(ax5);
    topoplot(topoLPP_mean.(condName), chanlocs, 'maplimits', clim, args{:});
    title(ax5, ttl)

end

xlabel(tA, 'Time (ms)');

% ============== Panel B layout: 3 rows (groups) x 3 cols (conditions) =================

tB = tiledlayout(pB, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for r = 1: 3

    grpName = cfg.msGroups{r};

    for c = 1: 3

        condName = cfg.condsLower{c};
        ax = nexttile(tB, (r - 1) * 3 + c);

        plot_ms_heatmap(ax, msProb.(grpName).(condName), cfg.msColors, {'A','B','C','D','E','F','G'}, cfg.plotMs, heatCLim);
        add_windows(ax, [cfg.winN200(1) cfg.winN200(2); cfg.winP300(1) cfg.winP300(2); cfg.winLPP(1) cfg.winLPP(2)], '');

        if r == 1
            title(ax, cfg.conds{c});
        end

        if c == 1
            ylabel(ax, grpName);
        end

        if r == 3
            xlabel(ax, 'Time (ms)');
        end
    end
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.TickDirection = 'out';
cb.Label.String = 'Label probability';

end

function cfg = set_defaults(cfg)

if ~isfield(cfg, 'conds')
    cfg.conds = {'Negative', 'Neutral', 'Positive'};
end

cfg.condsLower = {'negative', 'neutral', 'positive'};

if ~isfield(cfg, 'baselineMs')
    cfg.baselineMs = [-250 0];
end

if ~isfield(cfg, 'plotMs')
    cfg.plotMs = [50 1000];
end

if ~isfield(cfg, 'winN200')
    cfg.winN200 = [180 300];
end

if ~isfield(cfg, 'winP300')
    cfg.winP300 = [300 500];
end

if ~isfield(cfg, 'winLPP')
    cfg.winLPP = [500 1000];
end

if ~isfield(cfg, 'roiN200')
    cfg.roiN200 = {'Fz', 'FC1', 'FC2', 'FC3', 'FC4', 'Cz', 'AFF1h', 'AFF2h'};
end

if ~isfield(cfg, 'roiP300_LPP')
    cfg.roiP300_LPP = {'CPz', 'CP1', 'CP2', 'Pz', 'P1', 'P2', 'POz', 'PO3', 'PO4'};
end

if ~isfield(cfg, 'groupColors')
    cfg.groupColors = [0.02 0.41 1.00; 1 0 0; 0 0.60 0];
end

if ~isfield(cfg, 'msColors')
    cfg.msColors = lines(7);
end

cfg.msGroups = {'HC', 'Siblings', 'BD'};

end

function [chanlocs, timesAll] = get_chanlocs_times(France, refCond)

EEG = France.HC{1}.(refCond).EEG_Data;
chanlocs = EEG.chanlocs;
timesAll = EEG.times;

end

function idx = local_chan_idx(chanlocs, labels)

allLabs = string({chanlocs.labels});
idx = zeros(1, numel(labels));

for i = 1:numel(labels)

    hit = find(strcmpi(allLabs, labels{i}), 1, 'first');

    if isempty(hit)
        error('Channel label not found: %s', labels{i});
    end

    idx(i) = hit;
end

end

function out = compute_erp_by_group(groupCells, condName, roiIdx, baselineIdx, plotIdx)

out = struct();
out.sub = cell(1, numel(groupCells));

for g = 1: numel(groupCells)
    subjCell = groupCells{g};
    X = erp_matrix(subjCell, condName, roiIdx, baselineIdx, plotIdx);
    out.sub{g} = X;
end

end

function X = erp_matrix(subjCell, condName, roiIdx, baselineIdx, plotIdx)

nSub = numel(subjCell);
X = nan(nSub, sum(plotIdx));

for s = 1:nSub

    EEG = subjCell{s}.(condName).EEG_Data;
    data = double(EEG.data);

    base = mean(data(:, baselineIdx, :), 2);
    data = data - base;

    roi = mean(data(roiIdx, plotIdx, :), 1);
    roi = squeeze(roi);
    roi = mean(roi, 2);

    X(s, :) = roi.';
end

end

function topo = compute_topo_at_sample(subjCell, condName, baselineIdx, timeIdx)

nSub  = numel(subjCell);
nChan = size(subjCell{1}.(condName).EEG_Data.data, 1);
acc   = zeros(nChan, nSub);

for s = 1: nSub

    EEG  = subjCell{s}.(condName).EEG_Data;
    data = double(EEG.data);

    base = mean(data(:, baselineIdx, :), 2);
    data = data - base;

    v = squeeze(mean(data(:, timeIdx, :), 3));
    acc(:, s) = v(:);
end

topo = mean(acc, 2);

end

function topo = compute_topo_mean_window(subjCell, condName, baselineIdx, winIdx)

nSub  = numel(subjCell);
nChan = size(subjCell{1}.(condName).EEG_Data.data, 1);
acc   = zeros(nChan, nSub);

for s = 1: nSub

    EEG  = subjCell{s}.(condName).EEG_Data;
    data = double(EEG.data);

    base = mean(data(:, baselineIdx, :), 2);
    data = data - base;

    seg = mean(data(:, winIdx, :), 2);     % [chan x 1 x trials]
    seg = squeeze(seg);                    % [chan x trials]
    seg = mean(seg, 2);                    % [chan x 1] mean across trials

    acc(:, s) = seg(:);

end

topo = mean(acc, 2);

end

function clim = get_clim_peak_struct(topoPeakStruct, conds)

vals = [];

for c = 1: numel(conds)
    condName = conds{c};

    for g = 1: numel(topoPeakStruct.(condName))
        vals = [vals; topoPeakStruct.(condName){g}(:)];
    end

end

mx = max(abs(vals));
clim = [-mx mx];

end

function plot_erp(ax, t, erpStruct, groupColors)

hold(ax, 'on');
set(ax, 'Color', 'w');

for g = 1: numel(erpStruct.sub)

    X = erpStruct.sub{g};
    mu = mean(X, 1, 'omitnan');
    n = size(X, 1);

    se = std(X, 0, 1, 'omitnan') ./ sqrt(n);
    tc = tinv(0.975, n - 1);

    lo = mu - tc .* se;
    hi = mu + tc .* se;

    col = groupColors(g, :);

    fill(ax, [t fliplr(t)], [lo fliplr(hi)], col, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(ax, t, mu, '-', 'Color', col, 'LineWidth', 2, 'HandleVisibility', 'off');

end

xlim(ax, [t(1) t(end)]);
set(ax, 'XTick', [200 400 600 800 1000])
grid(ax, 'off');
set(ax, 'Box', 'off');

end

function add_windows(ax, wins, labelText)

yl = ylim(ax);

if isvector(wins)
    wins = reshape(wins, 1, 2);
end

for i = 1:size(wins, 1)
    xline(ax, wins(i, 1), ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    xline(ax, wins(i, 2), ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
end

ylim(ax, yl);

end

function msProb = compute_ms_prob(outputStats, msGroups, condsLower, plotMs)

tIdx = (plotMs(1): plotMs(2)) + 1;
tLen = numel(tIdx);

msProb = struct();

for g = 1: numel(msGroups)

    grpName = msGroups{g};

    for c = 1:numel(condsLower)

        condName = condsLower{c};

        subjArr = get_ms_subjects(outputStats, grpName, condName);
        P = zeros(7, tLen);
        nSub = numel(subjArr);

        for s = 1: nSub

            L = double(subjArr(s).MSClass);
            L = L(tIdx, :);

            for k = 1: 7
                P(k, :) = P(k, :) + mean(L == k, 2).';
            end

        end

        P = P ./ nSub;
        msProb.(grpName).(condName) = P;

    end
end

end

function subjArr = get_ms_subjects(outputStats, grpName, condName)

switch grpName

    case 'HC'
        subjArr = outputStats.HC.(condName);

    case 'Siblings'
        subjArr = outputStats.Siblings.(condName);

    case 'BD'
        a = outputStats.BP_I_Depressed.(condName);
        b = outputStats.BP_I_Euthymic.(condName);
        c = outputStats.BP_II_Depressed.(condName);
        d = outputStats.BP_II_Euthymic.(condName);
        subjArr = [a(:); b(:); c(:); d(:)];

    otherwise
        error('Unknown ms group: %s', grpName);

end

end

function plot_ms_heatmap(ax, P, msColors, mapLabels, plotMs, heatCLim)

t = plotMs(1): plotMs(2);
imagesc(ax, t, 1: 7, P);

set(ax, 'YDir', 'normal');
set(ax, 'Color', 'w');

colormap(ax, turbo);

if nargin < 6 || isempty(heatCLim)
    caxis(ax, [0 1]);
else
    caxis(ax, heatCLim);
end

xlim(ax, [t(1) t(end)]);
set(ax, 'XTick', [200 400 600 800 1000])
ylim(ax, [0.5 7.5]);
set(ax, 'Box', 'off');

if ~isempty(mapLabels)
    yticks(ax, 1: 7);
    yticklabels(ax, mapLabels);
else
    yticks(ax, 1: 7);
    yticklabels(ax, {'A', 'B', 'C', 'D', 'E', 'F', 'G'});
end

end