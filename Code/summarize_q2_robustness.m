function [overall, byERP] = summarize_q2_robustness(Base, Sens, ERPs, Conds, Groups)

% SUMMARIZE_Q2_ROBUSTNESS
% Compare sensitivity pipeline vs baseline for Q2 (target vs HC edge-wise rate ratios).
%
% Primary:
%   - Among BASELINE-significant edges: direction retention and significance retention.
%
% Added:
%   - Jaccard overlap, lost/gained counts
%   - Percent agreement and Cohen's kappa for edge significance
%   - Effect-size concordance on log(RR) (or log(IRR) fallback)
%
% IMPORTANT REVISION:
%   Jaccard and kappa are now summarized in TWO ways:
%     (1) "All cells" includes empty-empty cells (both pipelines have no sig edges)
%     (2) "Non-empty cells" excludes cells where union(sigB, sigS) is empty.
%   This avoids misleading medians of 1 driven by empty-empty cells.

alpha = 0.05;

edgeMask = ~eye(7);

total_base = 0;
retain_dir = 0;
retain_sig = 0;

cellStats = [];

per_erp = struct();

for e = 1: numel(ERPs)

    erp = ERPs{e};

    bE = 0;
    dE = 0;
    sE = 0;

    cellRows = [];

    for c = 1: numel(Conds)

        cond = Conds{c};

        panelB = safeget(Base, [], erp, cond);
        panelS = safeget(Sens, [], erp, cond);

        if isempty(panelB) || isempty(panelS)
            continue;
        end

        tgs = fieldnames(panelB);
        tgs = tgs(ismember(tgs, Groups) & ~strcmp(tgs, 'FDR_SCOPE'));

        for k = 1: numel(tgs)

            tgt = tgs{k};

            B = safeget(panelB, [], tgt);
            S = safeget(panelS, [], tgt);

            if isempty(B) || isempty(S)
                continue;
            end

            if ~isfield(B, 'sig') || ~isfield(S, 'sig')
                continue;
            end

            ratioFieldB = pick_ratio_field(B);
            ratioFieldS = pick_ratio_field(S);

            if strlength(ratioFieldB) == 0 || strlength(ratioFieldS) == 0
                continue;
            end

            RB = B.(ratioFieldB);
            RS = S.(ratioFieldS);

            sigB = logical(B.sig) & edgeMask;
            sigS = logical(S.sig) & edgeMask;

            % ---- Baseline significant edges: retain direction + retain significance ----

            if any(sigB(:))

                idx = find(sigB);

                for q = 1: numel(idx)

                    total_base = total_base + 1;
                    bE = bE + 1;

                    [i, j] = ind2sub([7 7], idx(q));

                    dirB = dir_from_ratio(RB(i, j));
                    dirS = dir_from_ratio(RS(i, j));

                    if ~isnan(dirB) && ~isnan(dirS) && dirB == dirS
                        retain_dir = retain_dir + 1;
                        dE = dE + 1;
                    end

                    if sigS(i, j)
                        retain_sig = retain_sig + 1;
                        sE = sE + 1;
                    end

                end

            end

            % ---- Added per-cell robustness stats ----

            [jacc, lostN, gainedN, agreePct, kappaVal, unionN] = sig_overlap_stats(sigB, sigS);
            [rho, rmseLog, madLog] = effect_concordance(RB, RS, edgeMask);

            cellRows = [cellRows; struct( ...
                'ERP', string(erp), ...
                'Condition', string(cond), ...
                'Target', string(tgt), ...
                'Jaccard', jacc, ...
                'UnionEdges', unionN, ...
                'LostEdges', lostN, ...
                'GainedEdges', gainedN, ...
                'AgreePct', agreePct, ...
                'Kappa', kappaVal, ...
                'RhoLog', rho, ...
                'RMSELog', rmseLog, ...
                'MADLog', madLog ...
                )];

        end
    end

    if isempty(cellRows)
        per_erp.(erp) = struct( ...
            'BaselineEdges', bE, ...
            'DirectionRetained', dE, ...
            'SignificanceRetained', sE, ...
            'DirPct', pct(dE, bE), ...
            'SigPct', pct(sE, bE), ...
            'CellsCompared', 0);
    else
        Tcell = struct2table(cellRows);
        cellStats = [cellStats; Tcell];

        % All-cells summaries (can be dominated by empty-empty cells)
        j_all = median(Tcell.Jaccard, 'omitnan');
        k_all = median(Tcell.Kappa, 'omitnan');

        % Non-empty cell summaries (recommended for interpretability)
        nonEmpty = (Tcell.UnionEdges > 0);
        if any(nonEmpty)
            j_ne = median(Tcell.Jaccard(nonEmpty), 'omitnan');
            k_ne = median(Tcell.Kappa(nonEmpty), 'omitnan');
            n_ne = sum(nonEmpty);
            n_kappa_def = sum(nonEmpty & isfinite(Tcell.Kappa));
        else
            j_ne = NaN;
            k_ne = NaN;
            n_ne = 0;
            n_kappa_def = 0;
        end

        per_erp.(erp) = struct( ...
            'BaselineEdges', bE, ...
            'DirectionRetained', dE, ...
            'SignificanceRetained', sE, ...
            'DirPct', pct(dE, bE), ...
            'SigPct', pct(sE, bE), ...
            'CellsCompared', height(Tcell), ...
            'Jaccard_Median_AllCells', j_all, ...
            'Kappa_Median_AllCells', k_all, ...
            'Jaccard_Median_NonEmpty', j_ne, ...
            'Kappa_Median_NonEmpty', k_ne, ...
            'NonEmptyCells', n_ne, ...
            'KappaDefined_NonEmptyCells', n_kappa_def, ...
            'RhoLog_Median', median(Tcell.RhoLog, 'omitnan'), ...
            'RMSELog_Median', median(Tcell.RMSELog, 'omitnan'), ...
            'MADLog_Median', median(Tcell.MADLog, 'omitnan'), ...
            'LostEdges_Total', sum(Tcell.LostEdges, 'omitnan'), ...
            'GainedEdges_Total', sum(Tcell.GainedEdges, 'omitnan') ...
            );

    end

end

overall = struct();
overall.BaselineEdges = total_base;
overall.DirectionRetained = retain_dir;
overall.SignificanceRetained = retain_sig;
overall.DirPct = pct(retain_dir, total_base);
overall.SigPct = pct(retain_sig, total_base);

if ~isempty(cellStats)

    overall.CellsCompared = height(cellStats);

    % All-cells summaries
    overall.Jaccard_Median_AllCells = median(cellStats.Jaccard, 'omitnan');
    overall.Kappa_Median_AllCells = median(cellStats.Kappa, 'omitnan');

    % Non-empty cells summaries
    nonEmpty = (cellStats.UnionEdges > 0);
    overall.NonEmptyCells = sum(nonEmpty);

    if any(nonEmpty)
        overall.Jaccard_Median_NonEmpty = median(cellStats.Jaccard(nonEmpty), 'omitnan');
        overall.Kappa_Median_NonEmpty = median(cellStats.Kappa(nonEmpty), 'omitnan');
        overall.KappaDefined_NonEmptyCells = sum(nonEmpty & isfinite(cellStats.Kappa));
    else
        overall.Jaccard_Median_NonEmpty = NaN;
        overall.Kappa_Median_NonEmpty = NaN;
        overall.KappaDefined_NonEmptyCells = 0;
    end

    overall.RhoLog_Median = median(cellStats.RhoLog, 'omitnan');
    overall.RMSELog_Median = median(cellStats.RMSELog, 'omitnan');
    overall.MADLog_Median = median(cellStats.MADLog, 'omitnan');
    overall.LostEdges_Total = sum(cellStats.LostEdges, 'omitnan');
    overall.GainedEdges_Total = sum(cellStats.GainedEdges, 'omitnan');
end

% ---- byERP table ----
% Add columns for all-cells and non-empty medians.

byERP = table('Size', [numel(ERPs) 14], ...
    'VariableTypes', {'string','double','double','double','string', ...
                      'double','double','double','double','double', ...
                      'double','double','double','double'}, ...
    'VariableNames', {'ERP','BaselineEdges','DirectionRetained','SignificanceRetained','Percentages', ...
                      'JaccardMed_All','KappaMed_All','JaccardMed_NonEmpty','KappaMed_NonEmpty','NonEmptyCells', ...
                      'RhoLogMed','RMSELogMed','MADLogMed','CellsCompared'});

for e = 1: numel(ERPs)

    erp = ERPs{e};
    S = per_erp.(erp);

    byERP.ERP(e) = string(erp);
    byERP.BaselineEdges(e) = S.BaselineEdges;
    byERP.DirectionRetained(e) = S.DirectionRetained;
    byERP.SignificanceRetained(e) = S.SignificanceRetained;
    byERP.Percentages(e) = sprintf('Dir=%s, Sig=%s', pct(S.DirectionRetained, S.BaselineEdges), pct(S.SignificanceRetained, S.BaselineEdges));

    byERP.CellsCompared(e) = S.CellsCompared;

    if isfield(S, 'Jaccard_Median_AllCells')
        byERP.JaccardMed_All(e) = S.Jaccard_Median_AllCells;
        byERP.KappaMed_All(e) = S.Kappa_Median_AllCells;
        byERP.JaccardMed_NonEmpty(e) = S.Jaccard_Median_NonEmpty;
        byERP.KappaMed_NonEmpty(e) = S.Kappa_Median_NonEmpty;
        byERP.NonEmptyCells(e) = S.NonEmptyCells;
        byERP.RhoLogMed(e) = S.RhoLog_Median;
        byERP.RMSELogMed(e) = S.RMSELog_Median;
        byERP.MADLogMed(e) = S.MADLog_Median;
    else
        byERP.JaccardMed_All(e) = NaN;
        byERP.KappaMed_All(e) = NaN;
        byERP.JaccardMed_NonEmpty(e) = NaN;
        byERP.KappaMed_NonEmpty(e) = NaN;
        byERP.NonEmptyCells(e) = 0;
        byERP.RhoLogMed(e) = NaN;
        byERP.RMSELogMed(e) = NaN;
        byERP.MADLogMed(e) = NaN;
    end

end

end


% ========================= HELPERS (local) =========================

function S = safeget(S0, defaultVal, varargin)

S = defaultVal;

try
    for k = 1: numel(varargin)
        S0 = S0.(varargin{k});
    end
    S = S0;
catch
end

end


function f = pick_ratio_field(S)

f = "";

if isfield(S, 'RR')
    f = "RR";
elseif isfield(S, 'IRR')
    f = "IRR";
end

end


function d = dir_from_ratio(r)

d = NaN;

if ~isfinite(r) || r <= 0
    return;
end

x = log(r);

if x > 0
    d = +1;
elseif x < 0
    d = -1;
else
    d = 0;
end

end


function s = pct(a, b)

if ~isfinite(a) || ~isfinite(b) || b <= 0
    s = 'NA';
    return;
end

s = sprintf('%.1f%%', 100 * a / b);

end


function [jacc, lostN, gainedN, agreePct, kappaVal, unionN] = sig_overlap_stats(sigB, sigS)

sigB = logical(sigB(:));
sigS = logical(sigS(:));

intN = sum(sigB & sigS);
unionN = sum(sigB | sigS);

if unionN == 0
    % Empty-empty cell; define jaccard as 1 but mark unionN=0 so it can be excluded.
    jacc = 1;
else
    jacc = intN / unionN;
end

lostN = sum(sigB & ~sigS);
gainedN = sum(~sigB & sigS);

agreePct = 100 * mean(sigB == sigS);

kappaVal = cohens_kappa(sigB, sigS);

end


function kappaVal = cohens_kappa(a, b)

a = logical(a(:));
b = logical(b(:));

po = mean(a == b);

pa1 = mean(a);
pb1 = mean(b);

pe = pa1 * pb1 + (1 - pa1) * (1 - pb1);

den = (1 - pe);

if den <= 0
    kappaVal = NaN;
else
    kappaVal = (po - pe) / den;
end

end


function [rho, rmseLog, madLog] = effect_concordance(RB, RS, edgeMask)

rho = NaN;
rmseLog = NaN;
madLog = NaN;

try
    vB = log(RB(edgeMask));
    vS = log(RS(edgeMask));
catch
    return;
end

good = isfinite(vB) & isfinite(vS);

if sum(good) < 3
    return;
end

vB = vB(good);
vS = vS(good);

try
    rho = corr(vB, vS, 'Type', 'Spearman');
catch
    rho = NaN;
end

rmseLog = sqrt(mean((vB - vS).^2, 'omitnan'));
madLog = median(abs(vB - vS), 'omitnan');

end
