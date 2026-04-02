function plot_emm_significant(OUT_measure, measure, ALL_Mean)

% PLOT_EMM_SIGNIFICANT
% Plot EMMs only when at least one pairwise comparison survives BOTH:
% 1) the ERP-wise (omnibus) BH-FDR gate, and
% 2) the within-family (cell) BH-FDR gate (pairwise pFDR < alpha).
%
% Usage:
% plot_emm_significant(Out_Duration_base, 'Duration', ALL_Mean);
% plot_emm_significant(Out_Occurrence_base, 'Num_occurrence', ALL_Mean);
% plot_emm_significant(Out_Coverage_base, 'Coverage', ALL_Mean);

alpha = 0.05;

if isfield(OUT_measure, 'info') && isfield(OUT_measure.info, 'alphaFDR') && isfinite(OUT_measure.info.alphaFDR)
    alpha = OUT_measure.info.alphaFDR;
end

ERPs = OUT_measure.info.ERPs;
Micros = OUT_measure.info.Microstates;

% Resolve measure field (supports calling with 'Occurrence' or 'Num_occurrence').
measureField = measure;

if strcmpi(measure, 'Occurrence') && isfield(OUT_measure, 'Num_occurrence')
    measureField = 'Num_occurrence';
end

for e = 1: numel(ERPs)
    for m = 1: numel(Micros)

        R = OUT_measure.(measureField){e, m};

        if isempty(R)
            continue;
        end

        erp = ERPs{e};
        ms = Micros{m};

        % ---------- Interaction significant (ERP BH-FDR) ----------

        pI = safeget(R, 'tests', 'interaction', 'pFDR');

        if ~isempty(pI) && isfinite(pI) && pI < alpha && isfield(R, 'pairwise')

            % Pairwise: Group within Condition
            P = get_first_table(R.pairwise, {'Group_within_Condition', 'Group_by_Condition'});
            T = get_first_table(R.emm, {'Group_by_Condition'});

            if ~isempty(P) && ~isempty(T) && ismember('Condition', P.Properties.VariableNames) && ismember('pFDR', P.Properties.VariableNames)

                condLv = categories(T.Condition);

                for c = 1: numel(condLv)

                    thisCond = condLv{c};

                    maskT = T.Condition == thisCond;
                    maskP = string(P.Condition) == string(thisCond);

                    Ps = P(maskP, :);
                    sigMask = isfinite(Ps.pFDR) & (Ps.pFDR < alpha);

                    if any(sigMask)
                        ttl = sprintf('%s | %s × %s — Groups @ Condition=%s', measureField, erp, ms, thisCond);
                        make_barplot_with_pairs(T(maskT, :), 'Group', Ps(sigMask, :), 'Group_1', 'Group_2', ttl, 'group', ALL_Mean, alpha);
                    end

                end

            end

            % Pairwise: Condition within Group
            P = get_first_table(R.pairwise, {'Condition_within_Group', 'Condition_by_Group'});
            T = get_first_table(R.emm, {'Condition_by_Group'});

            if ~isempty(P) && ~isempty(T) && ismember('Group', P.Properties.VariableNames) && ismember('pFDR', P.Properties.VariableNames)

                grpLv = categories(T.Group);

                for g = 1: numel(grpLv)

                    thisGrp = grpLv{g};

                    maskT = T.Group == thisGrp;
                    maskP = string(P.Group) == string(thisGrp);

                    Ps = P(maskP, :);
                    sigMask = isfinite(Ps.pFDR) & (Ps.pFDR < alpha);

                    if any(sigMask)
                        ttl = sprintf('%s | %s × %s — Conditions @ Group=%s', measureField, erp, ms, thisGrp);
                        make_barplot_with_pairs(T(maskT, :), 'Condition', Ps(sigMask, :), 'Condition_1', 'Condition_2', ttl, 'condition', ALL_Mean, alpha);
                    end

                end

            end

            continue; % Skip main-effects plotting when interaction is significant.

        end

        % ---------- No interaction: main effects (ERP BH-FDR) and cell BH-FDR pairs ----------

        pG = safeget(R, 'tests', 'group', 'pFDR');

        if ~isempty(pG) && isfinite(pG) && pG < alpha && isfield(R, 'pairwise') && isfield(R.pairwise, 'Group')

            T = safeget(R, 'emm', 'Group');
            P = R.pairwise.Group;

            if ~isempty(T) && istable(P) && ismember('pFDR', P.Properties.VariableNames)

                sigMask = isfinite(P.pFDR) & (P.pFDR < alpha);

                if any(sigMask)
                    ttl = sprintf('%s | %s × %s — Groups (marginal)', measureField, erp, ms);
                    make_barplot_with_pairs(T, 'Group', P(sigMask, :), 'Group_1', 'Group_2', ttl, 'group', ALL_Mean, alpha);
                end

            end

        end

        pC = safeget(R, 'tests', 'condition', 'pFDR');

        if ~isempty(pC) && isfinite(pC) && pC < alpha && isfield(R, 'pairwise') && isfield(R.pairwise, 'Condition')

            T = safeget(R, 'emm', 'Condition');
            P = R.pairwise.Condition;

            if ~isempty(T) && istable(P) && ismember('pFDR', P.Properties.VariableNames)

                sigMask = isfinite(P.pFDR) & (P.pFDR < alpha);

                if any(sigMask)
                    ttl = sprintf('%s | %s × %s — Conditions (marginal)', measureField, erp, ms);
                    make_barplot_with_pairs(T, 'Condition', P(sigMask, :), 'Condition_1', 'Condition_2', ttl, 'condition', ALL_Mean, alpha);
                end

            end

        end

    end
end

end


function make_barplot_with_pairs(T, factorCol, P, leftCol, rightCol, ttl, modeTag, ALL_Mean, alpha)

% MAKE_BARPLOT_WITH_PAIRS
% Bar plot of back-transformed EMMs with 95% CIs, raw-data violin overlay (if available),
% and pairwise-significance brackets for the (already filtered) significant pairs.

gr_name = {'BP_I_Depressed', 'BP_I_Euthymic', 'BP_II_Depressed', 'BP_II_Euthymic', 'HC', 'Siblings'};
cond_name = {'Neutral', 'Negative', 'Positive'};

[measureVar, ~, ERP, micro] = parse_plot_title(ttl);

% -------- Build raw data for violins (context-aware) --------

[RawData, RawName] = build_raw_for_violin(ALL_Mean, ERP, micro, measureVar, factorCol, T);

% -------- Extract EMMs --------

labels = string(T.(factorCol));
y = T.Estimate_BT;
yL = T.CI_lo_BT;
yH = T.CI_hi_BT;

if ismember('Units', T.Properties.VariableNames)
    units = string(T.Units(1));
else
    units = "";
end

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

x = 1: 4: numel(labels) * 4;
lo = max(0, y - yL);
hi = max(0, yH - y);

fig = figure('Color', 'w', 'Position', [100 100 950 520]);
b = bar(x, y, 0.4);
b.FaceColor = 0.75 * [1 1 1];
b.EdgeColor = 0.2 * [0 0 0];

hold on;

errorbar(x, y, lo, hi, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

set(gca, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 30, 'FontSize', 11, 'Box', 'off');
ylabel(sprintf('EMM (%s)', units));
title(ttl, 'Interpreter', 'none');

% -------- Violin overlay (optional; will not error if function missing) --------

try
    if ~isempty(RawData) && ~isempty(RawName)

        if strcmp(factorCol, 'Group')
            daviolinplot_microstate_task(RawData, 'groups', RawName, 'boxwidth', 2, 'colors', repmat(0.6 * [1 1 1], 6, 1), 'box', 0, 'outliers', 0, 'violinwidth', 3.5);
            set(gca, 'XLim', [-0.5 23.5]);
        elseif strcmp(factorCol, 'Condition')
            daviolinplot_microstate_task(RawData, 'groups', RawName, 'boxwidth', 2, 'colors', repmat(0.6 * [1 1 1], 3, 1), 'box', 0, 'outliers', 0, 'violinwidth', 2);
            set(gca, 'XLim', [-0.5 12.5]);
        end

    end
catch
    % If the violin function is not available, silently continue.
end

% -------- Pairwise brackets for significant pairs only --------
% P should already contain only pFDR < alpha rows, but we enforce it again for safety.

if istable(P) && ismember('pFDR', P.Properties.VariableNames)
    sigMask = isfinite(P.pFDR) & (P.pFDR < alpha);
    P = P(sigMask, :);
end

if ~isempty(P) && istable(P)

    [P, leftColResolved, rightColResolved] = resolve_pair_cols(P, factorCol);

    if ismember(leftColResolved, P.Properties.VariableNames) && ismember(rightColResolved, P.Properties.VariableNames)
        annotate_pairs(P, leftColResolved, rightColResolved, x, labels, yH, alpha);
    end

end

end


function annotate_pairs(P, leftCol, rightCol, x, labels, yH, alpha)

% Draw simple brackets with pFDR labels for each row of P.

lab = string(labels(:));
xMap = containers.Map();

for i = 1: numel(lab)
    xMap(char(lab(i))) = x(i);
end

yTop = max(yH);

yr = max(yTop - min(0, min(yH)), 1);
base = yTop + 0.05 * yr;
step = 0.06 * yr;

k = 0;

for r = 1: height(P)

    L1 = string(P.(leftCol)(r));
    L2 = string(P.(rightCol)(r));

    if ~isKey(xMap, char(L1)) || ~isKey(xMap, char(L2))
        continue;
    end

    if ismember('pFDR', P.Properties.VariableNames)
        ptxt = fmt_p(P.pFDR(r));
    else
        ptxt = "";
    end

    x1 = xMap(char(L1));
    x2 = xMap(char(L2));

    if x1 == x2
        continue;
    end

    k = k + 1;
    y = base + (k - 1) * step;

    plot([x1 x1 x2 x2], [y y + 0.02 * yr y + 0.02 * yr y], 'k', 'LineWidth', 1.2);
    text(mean([x1 x2]), y + 0.025 * yr, sprintf('pFDR %s', ptxt), 'HorizontalAlignment', 'center', 'FontSize', 10);

end

end


function [RawData, RawName] = build_raw_for_violin(ALL_Mean, ERP, micro, measureVar, factorCol, T)

% Build raw vectors for violin overlay. Uses context if T contains the other factor.
% - Group within Condition: use that Condition only.
% - Condition within Group: use that Group only.
% - Marginal Group: average across conditions within subject.
% - Marginal Condition: pool subjects across groups, within each condition.

RawData = [];
RawName = [];

S = ALL_Mean.(ERP).(micro);

subj = S.Number(:);
grp = string(S.Group(:));
cond = string(S.Condition(:));
val = S.(measureVar)(:);

Tb = table(subj, grp, cond, val, 'VariableNames', {'Subject', 'Group', 'Condition', 'Value'});

groupLevels = ["BP_I_Depressed", "BP_I_Euthymic", "BP_II_Depressed", "BP_II_Euthymic", "HC", "Siblings"];
condLevels = ["Neutral", "Negative", "Positive"];

if strcmpi(factorCol, 'Group')

    useCond = "";

    if ismember('Condition', T.Properties.VariableNames)
        u = unique(string(T.Condition));
        if isscalar(u)
            useCond = u(1);
        end
    end

    allVals = [];
    allNames = [];

    for i = 1: numel(groupLevels)

        g = groupLevels(i);
        Tg = Tb(Tb.Group == g, :);

        if strlength(useCond) > 0
            Tg = Tg(Tg.Condition == useCond, :);
            v = Tg.Value;
        else
            % Average within subject across conditions (subject-level mean).
            v = groupsummary(Tg, 'Subject', 'mean', 'Value');
            v = v.mean_Value;
        end

        allVals = [allVals; v(:)];
        allNames = [allNames; repmat(g, numel(v), 1)];

    end

    RawData = allVals;
    RawName = allNames;

elseif strcmpi(factorCol, 'Condition')

    useGroup = "";

    if ismember('Group', T.Properties.VariableNames)
        u = unique(string(T.Group));
        if isscalar(u)
            useGroup = u(1);
        end
    end

    allVals = [];
    allNames = [];

    for i = 1: numel(condLevels)

        c = condLevels(i);

        Tc = Tb(Tb.Condition == c, :);

        if strlength(useGroup) > 0
            Tc = Tc(Tc.Group == useGroup, :);
        end

        v = Tc.Value;

        allVals = [allVals; v(:)];
        allNames = [allNames; repmat(c, numel(v), 1)];

    end

    RawData = allVals;
    RawName = allNames;

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

% Returns [] if any field is missing.

val = [];

try
    for k = 1: numel(varargin)
        S = S.(varargin{k});
    end
    val = S;
catch
    val = [];
end

end


function T = get_first_table(S, candidates)

T = [];

if isempty(S) || ~isstruct(S)
    return;
end

for k = 1: numel(candidates)
    fn = candidates{k};
    if isfield(S, fn) && istable(S.(fn)) && ~isempty(S.(fn))
        T = S.(fn);
        return;
    end
end

end


function [P, leftCol, rightCol] = resolve_pair_cols(P, factorCol)

c = P.Properties.VariableNames;

if all(ismember({'Level1', 'Level2'}, c))
    leftCol = 'Level1';
    rightCol = 'Level2';
    return;
end

if strcmpi(factorCol, 'Group')

    if all(ismember({'Group_1', 'Group_2'}, c))
        leftCol = 'Group_1';
        rightCol = 'Group_2';
        return;
    end

elseif strcmpi(factorCol, 'Condition')

    if all(ismember({'Condition_1', 'Condition_2'}, c))
        leftCol = 'Condition_1';
        rightCol = 'Condition_2';
        return;
    end

end

candLeft = c(endsWith(c, '_1'));
candRight = c(endsWith(c, '_2'));

if ~isempty(candLeft) && ~isempty(candRight)
    leftCol = candLeft{1};
    rightCol = candRight{1};
    return;
end

% New format: a single contrast label column (e.g., P.Group = "A - B").
% Parse it into Level1 / Level2.

if ismember(factorCol, c)

    lab = string(P.(factorCol));
    L1 = strings(height(P), 1);
    L2 = strings(height(P), 1);

    for r = 1: height(P)

        s = strtrim(lab(r));

        if strlength(s) == 0
            continue;
        end

        parts = split(s, " - ");

        if numel(parts) < 2
            parts = regexp(s, '\s+vs\s+', 'split');
        end

        if numel(parts) >= 2
            L1(r) = strtrim(parts{1});
            L2(r) = strtrim(parts{2});
        end

    end

    if any(strlength(L1) > 0) && any(strlength(L2) > 0)

        P.Level1 = L1;
        P.Level2 = L2;
        leftCol = 'Level1';
        rightCol = 'Level2';
        return;
    end

end

error('Could not resolve pairwise label columns in the pairwise table.');

end


function [measureVar, measureLabel, ERP, microstate] = parse_plot_title(ttl)

if ~ischar(ttl) && ~isstring(ttl)
    error('ttl must be char or string');
end

ttl = string(ttl);

% 1) Remove anything after the em-dash (context like "— Groups (marginal)").
partsDash = regexp(ttl, '\s+—\s+', 'split');
head = strtrim(partsDash{1});

% 2) Split "<Measure> | <ERP> × <Microstate>"
partsPipe = regexp(head, '\s+\|\s+', 'split');

if numel(partsPipe) ~= 2
    error('Unexpected title format (pipe section). Title: %s', ttl);
end

measureLabel = strtrim(partsPipe{1});
rightSide = strtrim(partsPipe{2});

% 3) Split "<ERP> × <Microstate>" (accept unicode × or plain x)
partsTimes = regexp(rightSide, '\s+[×x]\s+', 'split', 'once');

if numel(partsTimes) ~= 2
    error('Unexpected title format (times section). Title: %s', ttl);
end

ERP = strtrim(partsTimes{1});
microstate = strtrim(partsTimes{2});

% 4) Map measure label to the data variable name
ml = lower(measureLabel);

switch ml
    case 'occurrence'
        measureVar = 'Num_occurrence';
    case 'num_occurrence'
        measureVar = 'Num_occurrence';
    case 'coverage'
        measureVar = 'Coverage';
    case 'duration'
        measureVar = 'Duration';
    otherwise
        % If you pass 'Num_occurrence' as the measure label in the title, use it directly.
        measureVar = measureLabel;
end

end

function s = fmt_p(p)

if isempty(p) || ~isfinite(p)
    s = 'NA';
elseif p < 0.001
    s = '<0.001';
else
    s = sprintf('%.3f', p);
end

end
