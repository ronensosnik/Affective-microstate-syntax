function plot_emm_significant(OUT_measure, measure, alphaFDR, ALL_Mean)

%% plot_emm_significant.m
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

if nargin < 3
    alphaFDR = 0.05;
end

ERPs  = OUT_measure.info.ERPs;
Micros = OUT_measure.info.Microstates;

for e = 1: numel(ERPs)
    for m = 1: numel(Micros)

        R = OUT_measure.(measure){e, m};

        if isempty(R)
            continue;
        end

        erp = ERPs{e};
        ms = Micros{m};

        % ---------- Interaction significant ----------

        pI = safeget(R, 'tests', 'interaction', 'pFDR');

        if ~isnan(pI) && pI < alphaFDR && isfield(R, 'pairwise')

            % Groups within each Condition

            if isfield(R.pairwise, 'Group_by_Condition') && isfield(R, 'emm') && isfield(R.emm, 'Group_by_Condition')

                P = R.pairwise.Group_by_Condition;
                T = R.emm.Group_by_Condition;

                if ismember('Condition', P.Properties.VariableNames)

                    condLv = categories(T.Condition);

                    for c = 1: numel(condLv)
                        maskT = T.Condition == condLv{c};
                        maskP = P.Condition == condLv{c};
                        Ps = P(maskP, :);

                        if any(~isnan(Ps.pFDR) & Ps.pFDR < alphaFDR)
                            ttl = sprintf('%s | %s × %s — Groups @ Condition=%s', measure, erp, ms, condLv{c});
                            make_barplot_with_pairs(T(maskT, :), 'Group', Ps, 'Group_1', 'Group_2', ttl, 'group');
                        end
                    end
                end
            end

            if isfield(R.pairwise, 'Condition_by_Group') && isfield(R, 'emm') && isfield(R.emm, 'Condition_by_Group')

                P = R.pairwise.Condition_by_Group;
                T = R.emm.Condition_by_Group;

                if ismember('Group', P.Properties.VariableNames)

                    grpLv = categories(T.Group);

                    for g = 1: numel(grpLv)

                        maskT = T.Group == grpLv{g};
                        maskP = P.Group == grpLv{g};
                        Ps = P(maskP, :);

                        if any(~isnan(Ps.pFDR) & Ps.pFDR < alphaFDR)
                            ttl = sprintf('%s | %s × %s — Conditions @ Group=%s', measure, erp, ms, grpLv{g});
                            make_barplot_with_pairs(T(maskT, :), 'Condition', Ps, 'Condition_1', 'Condition_2', ttl, 'condition');
                        end
                    end
                end
            end

            continue; % skip to next (we don't also plot main effects)
        end

        % ---------- No interaction: significant main effects only ----------

        pG = safeget(R, 'tests', 'group', 'pFDR');

        if ~isnan(pG) && pG < alphaFDR && isfield(R, 'pairwise') && isfield(R.pairwise, 'Group')

            T = safeget(R, 'emm', 'Group');
            P = R.pairwise.Group;

            if ~isempty(T) && any(~isnan(P.pFDR) & P.pFDR < alphaFDR)
                ttl = sprintf('%s | %s × %s — Groups (marginal)', measure, erp, ms);
                make_barplot_with_pairs(T, 'Group', P, 'Group_1','Group_2', ttl, 'group', ALL_Mean);
            end
        end

        pC = safeget(R, 'tests', 'condition', 'pFDR');

        if ~isnan(pC) && pC < alphaFDR && isfield(R,'pairwise') && isfield(R.pairwise, 'Condition')
            T = safeget(R, 'emm', 'Condition');
            P = R.pairwise.Condition;

            if ~isempty(T) && any(~isnan(P.pFDR) & P.pFDR < alphaFDR)
                ttl = sprintf('%s | %s × %s — Conditions (marginal)', measure, erp, ms);
                make_barplot_with_pairs(T, 'Condition', P, 'Condition_1', 'Condition_2', ttl, 'condition', ALL_Mean);
            end
        end

    end
end

end

function make_barplot_with_pairs(T, factorCol, P, leftCol, rightCol, ttl, modeTag, ALL_Mean)

gr_name = {'BP_I_Depressed', 'BP_I_Euthymic', 'BP_II_Depressed', 'BP_II_Euthymic', 'HC', 'Siblings'};
cond_name = {'Neutral', 'Negative', 'Positive'};

[measureVar, measureLabel, ERP, micro] = parse_plot_title(ttl);

if strcmp(factorCol, 'Group')

    for pop = 1: 6

        clear aa; aa = find(strcmp(ALL_Mean.(ERP).(micro).Group, gr_name{pop}))';

        clear bb; bb = find(strcmp(ALL_Mean.(ERP).(micro).Condition, 'Neutral'))';
        clear cc; cc  = find(strcmp(ALL_Mean.(ERP).(micro).Condition, 'Negative'))';
        clear dd; dd = find(strcmp(ALL_Mean.(ERP).(micro).Condition, 'Positive'))';

        Raw.data.(gr_name{pop}) = mean([ALL_Mean.(ERP).(micro).(measureVar)(intersect(aa, bb, 'stable')); ALL_Mean.(ERP).(micro).(measureVar)(intersect(aa, cc, 'stable')); ALL_Mean.(ERP).(micro).(measureVar)(intersect(aa, dd, 'stable'))]', 2);
        n = numel(Raw.data.(gr_name{pop}));

        for popp = 1: length(intersect(aa, bb, 'stable'))
            Raw.name.(gr_name{pop}){popp} = gr_name{pop};
        end

    end

    clear cells; cells = cellfun(@(g) Raw.data.(g)(:), gr_name, 'UniformOutput', false);
    RawData = vertcat(cells{:});

    clear cells; cells = cellfun(@(g) Raw.name.(g)(:), gr_name, 'UniformOutput', false);
    RawName = vertcat(cells{:});

elseif strcmp(factorCol, 'Condition')

    for pop = 1: 3

        clear aa; aa = find(strcmp(ALL_Mean.(ERP).(micro).Condition, cond_name{pop}))';

        clear bb; bb = find(strcmp(ALL_Mean.(ERP).(micro).Group, 'BP_I_Depressed'))';
        clear cc; cc  = find(strcmp(ALL_Mean.(ERP).(micro).Group, 'BP_I_Euthymic'))';
        clear dd; dd = find(strcmp(ALL_Mean.(ERP).(micro).Group, 'BP_II_Depressed'))';
        clear ee; ee  = find(strcmp(ALL_Mean.(ERP).(micro).Group, 'BP_II_Euthymic'))';
        clear ff;   ff    = find(strcmp(ALL_Mean.(ERP).(micro).Group, 'HC'))';
        clear gg; gg  = find(strcmp(ALL_Mean.(ERP).(micro).Group, 'Siblings'))';

        Raw.data.(cond_name{pop}) =  [ALL_Mean.(ERP).(micro).(measureVar)(intersect(aa, bb, 'stable'))'; ALL_Mean.(ERP).(micro).(measureVar)(intersect(aa, cc, 'stable'))'; ALL_Mean.(ERP).(micro).(measureVar)(intersect(aa, dd, 'stable'))'; ALL_Mean.(ERP).(micro).(measureVar)(intersect(aa, ee, 'stable'))'; ALL_Mean.(ERP).(micro).(measureVar)(intersect(aa, ff, 'stable'))'; ALL_Mean.(ERP).(micro).(measureVar)(intersect(aa, gg, 'stable'))'];
        n = numel(Raw.data.(cond_name{pop}));

        for popp = 1: n
            Raw.name.(cond_name{pop}){popp} = cond_name{pop};
        end

    end

    clear cells; cells = cellfun(@(g) Raw.data.(g)(:), cond_name, 'UniformOutput', false);
    RawData = vertcat(cells{:});

    clear cells; cells = cellfun(@(g) Raw.name.(g)(:), cond_name, 'UniformOutput', false);
    RawName = vertcat(cells{:});

end

labels = string(T.(factorCol));
y        = T.EMM_BT;
yL      = T.EMM_BT_Low;
yH     = T.EMM_BT_High;
units  = string(T.Units(1));

if strcmpi(modeTag, 'group')

    wanted = ["BP_I_Depressed", "BP_I_Euthymic", "BP_II_Depressed", "BP_II_Euthymic", "HC", "Siblings"];
else
    wanted = ["Neutral", "Negative", "Positive"];
end

[ord, ~] = reorder(labels, wanted);
labels = labels(ord);
y = y(ord);
yL = yL(ord);
yH = yH(ord);

x  = 1: 4: numel(labels) * 4;
lo = max(0, y - yL);
hi = max(0, yH - y);

fig = figure('Color', 'w', 'Position', [100 100 950 520]);
b   = bar(x, y, 0.4);
b.FaceColor = 0.75 * [1 1 1];
b.EdgeColor = 0.2 * [0 0 0];

hold on;

errorbar(x, y, lo, hi, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 30, 'FontSize', 11, 'Box', 'off');
ylabel(sprintf('EMM (%s)', units));
title(ttl, 'Interpreter', 'none');

maxH = NaN;
drew = false;
[leftCol, rightCol] = resolve_pair_cols(P, factorCol);

if strcmp(factorCol, 'Group')
    h = daviolinplot_microstate_task(RawData, 'groups', RawName, 'boxwidth', 2, 'colors', [0.6 0.6 0.6; 0.6 0.6 0.6; 0.6 0.6 0.6; 0.6 0.6 0.6; 0.6 0.6 0.6; 0.6 0.6 0.6], 'box', 0, 'outliers', 0, 'violinwidth', 3.5);
    set(gca, 'XLim', [-0.5 23.5]);
elseif strcmp(factorCol, 'Condition')
    h = daviolinplot_microstate_task(RawData, 'groups', RawName, 'boxwidth', 2, 'colors', [0.6 0.6 0.6; 0.6 0.6 0.6; 0.6 0.6 0.6], 'box', 0, 'outliers', 0, 'violinwidth', 2);
    set(gca, 'XLim', [-0.5 12.5]);
end

end

function [ord, keepMask] = reorder(labels, wanted)

labels = string(labels);
wanted = string(wanted);
keepMask = ismember(wanted, labels);
present = wanted(keepMask);
ord = zeros(1, numel(present));

for i = 1: numel(present)
    ord(i) = find(labels == present(i), 1, 'first');
end
end

function RGB = col(c)

if isnumeric(c) && numel(c) == 3
    RGB = c(:)';
    return;
end

switch char(c)
    case 'm'
        RGB = [1 0 1];
    case 'b'
        RGB = [0 0 1];
    case 'g'
        RGB = [0 1 0];
    case 'k'
        RGB = [0 0 0];
    case 'r'
        RGB = [1 0 0];
    case 'c'
        RGB = [0 1 1];
    case 'y'
        RGB = [1 1 0];
    otherwise
        RGB = 0.75 * [1 1 1];
end
end

function val = safeget(S, varargin)

val = NaN;

try
    for k=1: numel(varargin)
        S = S.(varargin{k});
    end
    val = S;
catch
end
end

function [leftCol, rightCol] = resolve_pair_cols(P, factorCol)

c = P.Properties.VariableNames;

if all(ismember({'Level1', 'Level2'}, c))
    leftCol   = 'Level1';
    rightCol = 'Level2';
    return;
end

if strcmpi(factorCol, 'Group')
    if all(ismember({'Group_1', 'Group_2'}, c))
        leftCol   = 'Group_1';
        rightCol = 'Group_2';
        return;
    end

elseif strcmpi(factorCol, 'Condition')
    if all(ismember({'Condition_1', 'Condition_2'}, c))
        leftCol   = 'Condition_1';
        rightCol = 'Condition_2';
        return;
    end
end

candLeft  = c(endsWith(c, '_1'));
candRight = c(endsWith(c, '_2'));

if ~isempty(candLeft) && ~isempty(candRight)
    leftCol  = candLeft{1};
    rightCol = candRight{1};
    return;
end

error('Could not resolve pairwise label columns in the pairwise table.');
end

function [measureVar, measureLabel, ERP, microstate] = parse_plot_title(ttl)

if ~ischar(ttl) && ~isstring(ttl)
    error('ttl must be char or string');
end

ttl = string(ttl);

% 1) Remove anything after the em-dash (context like "— Groups (marginal)")

partsDash = regexp(ttl, '\s+—\s+', 'split');   % split on em-dash
head = strtrim(partsDash{1});                  % keep left part

% 2) Split "<Measure> | <ERP> × <Microstate>"

partsPipe = regexp(head, '\s+\|\s+', 'split'); % split on pipe

if numel(partsPipe) ~= 2
    error('Unexpected title format (pipe section). Title: %s', ttl);
end
measureLabel = strtrim(partsPipe{1});
rightSide = strtrim(partsPipe{2});

% 3) Split "<ERP> × <Microstate>" (accept the unicode × or a plain x)

partsTimes = regexp(rightSide, '\s+[×x]\s+', 'split', 'once');

if numel(partsTimes) ~= 2
    error('Unexpected title format (times section). Title: %s', ttl);
end

ERP = strtrim(partsTimes{1});
microstate = strtrim(partsTimes{2});

% 4) Map measure label to the data variable name you want

switch lower(measureLabel)
    case 'occurrence'
        measureVar = 'Num_occurrence';
    case 'coverage'
        measureVar = 'Coverage';
    case 'duration'
        measureVar = 'Duration';
    otherwise
        measureVar = measureLabel;
end
end
