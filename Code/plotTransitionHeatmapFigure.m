function plotTransitionHeatmapFigure(TransStats, CompStats, ERPs, Conditions, Groups, varargin)

% plotTransitionHeatmapFigure: Matrix-view replacement for circular transition plots (Q1 + Q2).
%
% Produces one figure with panels for ERPs:
% - Panel A: N200
% - Panel B: LPP
%
% Each ERP panel contains rows for groups with any significant result in that ERP.
% Within each row, there are 3 condition blocks (Negative / Neutral / Positive).
% Within each condition block:
% - Left heatmap: Q1 (within-group deviation from independence), shown as log(IRR) with FDR mask.
% - Right heatmap: Q2 (group vs HC), shown as log(RR) with FDR mask (NA for HC).
%
% Notes:
% - Row index = source microstate; column index = target microstate.
% - Only FDR-significant edges are shown (others are masked).
% - Color limits are symmetric about 0 and shared within each ERP panel.

p = inputParser;
addParameter(p, 'MicrostateLabels', {'A', 'B', 'C', 'D', 'E', 'F', 'G'});
addParameter(p, 'PanelOrder', {'N200', 'LPP', 'P300'});
addParameter(p, 'AlwaysInclude', {'N200', 'LPP'});
addParameter(p, 'AlphaFDR', 0.05);
parse(p, varargin{:});
opt = p.Results;

MicrostateLabels = opt.MicrostateLabels;

% Decide which ERP panels to show

wantERPs = opt.PanelOrder;
present = intersect(wantERPs, ERPs, 'stable');

panelERPs = {};

for i = 1: numel(present)
    panelERPs{end + 1} = present{i};
end

% Always include N200 and LPP if present

for i = 1: numel(opt.AlwaysInclude)

    erpNeed = opt.AlwaysInclude{i};

    if any(strcmp(panelERPs, erpNeed))
        continue;
    end

    if any(strcmp(present, erpNeed))
        panelERPs = [{erpNeed} panelERPs];
    end
    
end

panelERPs = unique(panelERPs, 'stable');

% Include P300 only if significant edges exist (unless not present)

if any(strcmp(panelERPs, 'P300'))
    if ~local_erp_has_any_sig(TransStats, CompStats, 'P300', Conditions, Groups, opt.AlphaFDR)
        panelERPs(strcmp(panelERPs, 'P300')) = [];
    end
end

% Keep only N200/LPP if requested ERPs are missing

if isempty(panelERPs)
    panelERPs = present;
end

nPanels = numel(panelERPs);

figure('Color', 'w', 'Units', 'normalized', 'Position', [0.05 0.05 0.90 0.88]);

% Panel layout (equal heights)

x0 = 0.05;
w0 = 0.90;
y0 = 0.05;
h0 = 0.90;
gap = 0.03;

if nPanels == 1
    panelPos = [x0 y0 w0 h0];
else
    hPanel = (h0 - gap * (nPanels - 1)) / nPanels;
    panelPos = zeros(nPanels, 4);

    for i = 1: nPanels
        panelPos(i, :) = [x0 y0 + (nPanels - i) * (hPanel + gap) w0 hPanel];
    end
end

panelLetters = {'A', 'B', 'C'};

for pIdx = 1: nPanels

    erp = panelERPs{pIdx};

    % Determine which groups get a row in this ERP panel

    groupsInPanel = local_groups_with_sig(TransStats, CompStats, erp, Conditions, Groups, opt.AlphaFDR);

    if isempty(groupsInPanel)
        groupsInPanel = {};
    end

    nRows = numel(groupsInPanel);

    % If no significant groups exist, show a single row with a note

    if nRows == 0
        nRows = 1;
        groupsInPanel = {'(no significant edges)'};
    end

    pan = uipanel('Units', 'normalized', 'Position', panelPos(pIdx, :), 'BorderType', 'none', 'BackgroundColor', 'w');

    annotation('textbox', [panelPos(pIdx, 1) - 0.045 panelPos(pIdx, 2) + panelPos(pIdx, 4) - 0.03 0.03 0.03], 'String', panelLetters{min(pIdx, numel(panelLetters))}, 'LineStyle', 'none', 'FontWeight', 'bold', 'FontSize', 14);

    titleBox = annotation('textbox', [panelPos(pIdx, 1) panelPos(pIdx, 2) + panelPos(pIdx, 4) - 0.03 panelPos(pIdx, 3) 0.03], 'String', erp, 'LineStyle', 'none', 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'left');
    titleBox.Color = 'w';

    t = tiledlayout(pan, nRows, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Compute symmetric color limits for this ERP panel (shared across all tiles in this panel)

    clim = local_panel_clim(TransStats, CompStats, erp, Conditions, groupsInPanel, Groups, opt.AlphaFDR);

    for r = 1: nRows

        grp = groupsInPanel{r};

        for c = 1: numel(Conditions)

            condName = Conditions{c};

            % -----------------------
            % Q1 heatmap (left)
            % -----------------------

            ax1 = nexttile(t, (r - 1) * 6 + (c - 1) * 2 + 1);

            if strcmp(grp, '(no significant edges)')
                axis(ax1, 'off');

                if r == 1 && c == 2
                    text(ax1, 0.5, 0.5, 'No FDR-significant edges', 'HorizontalAlignment', 'center');
                end
            else

                [M1, has1] = local_get_q1_matrix(TransStats, erp, grp, condName, opt.AlphaFDR);
                local_plot_matrix(ax1, M1, clim, MicrostateLabels);

                if r == 1
                    title(ax1, sprintf('%s Q1', condName), 'FontWeight', 'normal', 'FontSize', 9);
                end

                if c == 1
                    ylabel(ax1, local_pretty_group(grp), 'FontWeight', 'bold');
                end
            end

            % -------------------------
            % Q2 heatmap (right)
            % -------------------------

            ax2 = nexttile(t, (r - 1) * 6 + (c - 1) * 2 + 2);

            if strcmp(grp, '(no significant edges)')
                axis(ax2, 'off');
            else
                [M2, has2] = local_get_q2_matrix(CompStats, erp, grp, condName, opt.AlphaFDR);
                local_plot_matrix(ax2, M2, clim, MicrostateLabels);

                if r == 1
                    title(ax2, sprintf('%s Q2', condName), 'FontWeight', 'normal', 'FontSize', 9);
                end
            end

            % Tick label density control

            if ~strcmp(grp, '(no significant edges)')
                if r < nRows
                    ax1.XTickLabel = {};
                    ax2.XTickLabel = {};
                end

                if c > 1
                    ax1.YTickLabel = {};
                    ax2.YTickLabel = {};
                else
                    ax2.YTickLabel = {};
                end
            end

        end
    end

    try
        cb = colorbar(t, 'eastoutside');
    catch
        cb = colorbar('Location', 'eastoutside');

        % If available, dock the colorbar to the tiledlayout's east tile.

        try
            cb.Layout.Tile = 'east';
        catch
        end
    end

    cb.TickDirection = 'out';
    cb.Label.String = 'log ratio (Q1: log IRR, Q2: log RR)';
    cb.FontSize = 9;

end

end

function tf = local_erp_has_any_sig(TransStats, CompStats, erp, Conditions, Groups, AlphaFDR)

tf = false;

% Q1

if isfield(TransStats, erp)
    for g = 1: numel(Groups)
        grp = Groups{g};

        if ~isfield(TransStats.(erp), grp)
            continue;
        end

        for c = 1: numel(Conditions)
            condName = Conditions{c};

            if isfield(TransStats.(erp).(grp), condName)
                sig = TransStats.(erp).(grp).(condName).sig;

                if any(sig(:))
                    tf = true;
                    return;
                end
            end
        end
    end
end

% Q2

if isfield(CompStats, erp)
    for c = 1: numel(Conditions)

        condName = Conditions{c};

        if ~isfield(CompStats.(erp), condName)
            continue;
        end

        for g = 1: numel(Groups)

            grp = Groups{g};

            if strcmp(grp, 'HC')
                continue;
            end

            if isfield(CompStats.(erp).(condName), grp)
                sig = CompStats.(erp).(condName).(grp).sig;

                if any(sig(:))
                    tf = true;
                    return;
                end
            end
        end
    end
end

end

function groupsInPanel = local_groups_with_sig(TransStats, CompStats, erp, Conditions, Groups, AlphaFDR)

groupsInPanel = {};

for g = 1: numel(Groups)

    grp = Groups{g};
    hasSig = false;

    % Q1

    if isfield(TransStats, erp) && isfield(TransStats.(erp), grp)
        for c = 1: numel(Conditions)

            condName = Conditions{c};

            if isfield(TransStats.(erp).(grp), condName)
                sig = TransStats.(erp).(grp).(condName).sig;

                if any(sig(:))
                    hasSig = true;
                    break;
                end
            end
        end
    end

    % Q2

    if ~hasSig && ~strcmp(grp, 'HC') && isfield(CompStats, erp)
        for c = 1: numel(Conditions)

            condName = Conditions{c};

            if isfield(CompStats.(erp), condName) && isfield(CompStats.(erp).(condName), grp)

                sig = CompStats.(erp).(condName).(grp).sig;

                if any(sig(:))
                    hasSig = true;
                    break;
                end
            end
        end
    end

    if hasSig
        groupsInPanel{end + 1} = grp;
    end

end

end

function clim = local_panel_clim(TransStats, CompStats, erp, Conditions, groupsInPanel, GroupsAll, AlphaFDR)

vals = [];

for r = 1: numel(groupsInPanel)
    grp = groupsInPanel{r};

    for c = 1: numel(Conditions)
        condName = Conditions{c};

        [M1, has1] = local_get_q1_matrix(TransStats, erp, grp, condName, AlphaFDR);
        if has1
            vals = [vals; M1(isfinite(M1))];
        end

        [M2, has2] = local_get_q2_matrix(CompStats, erp, grp, condName, AlphaFDR);
        if has2
            vals = [vals; M2(isfinite(M2))];
        end
    end
end

if isempty(vals)
    mx = 1;
else
    mx = prctile(abs(vals), 70);

    if ~isfinite(mx) || mx <= eps
        mx = 1;
    end
end

clim = [-mx mx];

end

function [M, hasAny] = local_get_q1_matrix(TransStats, erp, grp, condName, AlphaFDR)

M = NaN(7);
hasAny = false;

if ~isfield(TransStats, erp)
    return;
end

if ~isfield(TransStats.(erp), grp)
    return;
end

if ~isfield(TransStats.(erp).(grp), condName)
    return;
end

IRR = TransStats.(erp).(grp).(condName).IRR;
sig = TransStats.(erp).(grp).(condName).sig;

if isempty(IRR) || isempty(sig)
    return;
end

M = log(IRR);
M(~sig) = NaN;
M(1: 8: end) = NaN;

hasAny = any(isfinite(M(:)));

end

function [M, hasAny] = local_get_q2_matrix(CompStats, erp, grp, condName, AlphaFDR)

M = NaN(7);
hasAny = false;

if strcmp(grp, 'HC')
    return;
end

if ~isfield(CompStats, erp)
    return;
end

if ~isfield(CompStats.(erp), condName)
    return;
end

if ~isfield(CompStats.(erp).(condName), grp)
    return;
end

IRR = CompStats.(erp).(condName).(grp).IRR;
sig = CompStats.(erp).(condName).(grp).sig;

if isempty(IRR) || isempty(sig)
    return;
end

M = log(IRR);
M(~sig) = NaN;
M(1: 8: end) = NaN;

hasAny = any(isfinite(M(:)));

end

function local_plot_matrix(ax, M, clim, MicrostateLabels)

hold(ax, 'on');
set(ax, 'Color', [1 1 1]);

h = imagesc(ax, M);
set(h, 'AlphaData', isfinite(M));

axis(ax, 'square');
xlim(ax, [0.5 7.5]);
ylim(ax, [0.5 7.5]);

set(ax, 'YDir', 'normal');

xticks(ax, 1: 7);
yticks(ax, 1: 7);

xticklabels(ax, MicrostateLabels);
yticklabels(ax, MicrostateLabels);

ax.TickLength = [0 0];
ax.FontSize = 7;

caxis(ax, clim);
colormap(ax, local_diverging_colormap(256));

end

function s = local_pretty_group(grp)

switch grp
    case 'BP_I_Depressed'
        s = 'BD-I Depressed';
    case 'BP_I_Euthymic'
        s = 'BD-I Euthymic';
    case 'BP_II_Depressed'
        s = 'BD-II Depressed';
    case 'BP_II_Euthymic'
        s = 'BD-II Euthymic';
    case 'HC'
        s = 'HC';
    case 'Siblings'
        s = 'Siblings';
    otherwise
        s = grp;
end

end

function cmap = local_diverging_colormap(n)

if nargin < 1
    n = 256;
end

n2 = floor(n / 2);

b = [linspace(0, 1, n2)' linspace(0, 1, n2)' ones(n2, 1)];
r = [ones(n - n2, 1) linspace(1, 0, n - n2)' linspace(1, 0, n - n2)'];

cmap = [b; r];

end