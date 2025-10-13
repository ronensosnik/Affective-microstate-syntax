function print_summary_measure(OUT_measure, alphaFDR)

%% print_summary_measure.m
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


% PRINT_SUMMARY_MEASURE
% Pretty, on-screen summary per ERP × Microstate for a single measure.
% Shows EMMs (back-transformed) and pairwise with directionality.
%
% Usage: print_summary_measure(Out_Coverage, 0.05)

if nargin < 2
    alphaFDR = 0.05;
end

ERPs     = OUT_measure.info.ERPs;
Micros   = OUT_measure.info.Microstates;
measure  = OUT_measure.info.measure;

barline  = repmat('-', 1, 78);
fprintf('\n%s\nSUMMARY (%s) — FDR alpha = %.3f\n%s\n', barline, measure, alphaFDR, barline);

for e = 1:numel(ERPs)
  for m = 1:numel(Micros)

    R = OUT_measure.(measure){e,m};

    if isempty(R)
        continue;
    end

    % Omnibus (already FDR-adjusted across the 21 tests)

    pI   = get_nested(R, {'tests', 'interaction', 'pFDR'});
    pG = get_nested(R, {'tests', 'group', 'pFDR'});
    pC = get_nested(R, {'tests', 'condition', 'pFDR'});

    % Model used + Age effect

    if strcmp(R.model_used, 'full_with_interaction')
        M = R.models.full_with_interaction;
    else
        if isfield(R.models,'no_interaction_REML')
            M = R.models.no_interaction_REML;
        else
            M = R.models.no_interaction;
        end
    end

    try
        A = anova(M, 'DFMethod', 'Satterthwaite');
        iAge = find(strcmp(A.Term,'Age_c'), 1);

        if ~isempty(iAge)
            FAge = A.FStat(iAge);
            pAge = A.pValue(iAge);
            df1 = A.DF1(iAge); 
            df2 = A.DF2(iAge);
        else
            FAge = NaN; 
            pAge = NaN; 
            df1 = NaN; 
            df2 = NaN;
        end

    catch
        FAge = NaN; 
        pAge = NaN;
        df1 = NaN;
        df2 = NaN;
    end

    header = sprintf('%s | %s × %s', measure, ERPs{e}, Micros{m});

    fprintf('\n%s\n%s\n', header, repmat('=', 1, numel(header)));
    fprintf('Model used: %s\n', strrep(R.model_used, '_',' '));
    fprintf('Age effect: F(%s,%s) = %s, p = %s  %s\n', fdf(df1), fdf(df2), fnum(FAge,3), fp(pAge), star(pAge < alphaFDR));
    fprintf('Omnibus (BH–FDR across 21): Interaction pFDR = %s %s | Group pFDR = %s %s | Condition pFDR = %s %s\n', fp(pI), star(pI<alphaFDR), fp(pG), star(pG<alphaFDR), fp(pC), star(pC<alphaFDR));

    % ================== Post-hoc reporting ==================

    if pI < alphaFDR

        % ----- Interaction present: simple-effects -----

        if isfield(R,'emm') && isfield(R.emm,'Group_by_Condition')
            T = R.emm.Group_by_Condition;
            fprintf('\nEMMs: Group within each Condition (%s)\n', units_for(measure));
            disp(T(:, [{'Group','Condition','EMM_BT','EMM_BT_Low','EMM_BT_High','Units'}]));
        end

        if isfield(R,'emm') && isfield(R.emm,'Condition_by_Group')
            T = R.emm.Condition_by_Group;
            fprintf('\nEMMs: Condition within each Group (%s)\n', units_for(measure));
            disp(T(:, [{'Condition','Group','EMM_BT','EMM_BT_Low','EMM_BT_High','Units'}]));
        end

        if isfield(R,'pairwise') && isfield(R.pairwise,'Group_within_Condition')
            fprintf('\nPairwise: Groups within each Condition (direction & CI)\n');
            print_pairs_simple(R.pairwise.Group_within_Condition, R.emm.Group_by_Condition, measure, 'Group','Condition', alphaFDR);
        end

        if isfield(R,'pairwise') && isfield(R.pairwise,'Condition_within_Group')
            fprintf('\nPairwise: Conditions within each Group (direction & CI)\n');
            print_pairs_simple(R.pairwise.Condition_within_Group, R.emm.Condition_by_Group, measure, 'Condition','Group', alphaFDR);
        end

    else

        % ----- No interaction: main effects if FDR-significant -----

        if pG < alphaFDR && isfield(R,'emm') && isfield(R.emm,'Group')
            fprintf('\nEMMs: Group (%s)\n', units_for(measure));
            disp(R.emm.Group(:, [{'Group','EMM_BT','EMM_BT_Low','EMM_BT_High','Units'}]));
        end

        if pC < alphaFDR && isfield(R,'emm') && isfield(R.emm,'Condition')
            fprintf('\nEMMs: Condition (%s)\n', units_for(measure));
            disp(R.emm.Condition(:, [{'Condition','EMM_BT','EMM_BT_Low','EMM_BT_High','Units'}]));
        end

        if isfield(R,'pairwise') && isfield(R.pairwise,'Group') && pG < alphaFDR
            fprintf('\nPairwise: Group (direction & CI)\n');
            print_pairs_main(R.pairwise.Group, R.emm.Group, measure, 'Group', alphaFDR);
        end

        if isfield(R, 'pairwise') && isfield(R.pairwise,'Condition') && pC < alphaFDR
            fprintf('\nPairwise: Condition (direction & CI)\n');
            print_pairs_main(R.pairwise.Condition, R.emm.Condition, measure, 'Condition', alphaFDR);
        end
    end
  end
end

fprintf('%s\n', barline);
end

% ---------------- local printers / utils ----------------

function print_pairs_main(P, T, measure, factorName, alphaFDR)

if ~ismember('pFDR', P.Properties.VariableNames)
    P.pFDR = bh_adjust(P.pValue);
end

sig = P.pFDR < alphaFDR & ~isnan(P.pFDR);

if ~any(sig)
    fprintf('  No pairwise survived FDR.\n'); 
    return;
end

u = units_for(measure);
rows = find(sig);

for r = rows(:).'
    L1 = string(P.Level1(r));
    L2 = string(P.Level2(r));
    E1 = pick_emm(T, factorName, L1);
    E2 = pick_emm(T, factorName, L2);

    switch lower(measure)
        case 'duration'
            d  = P.Delta(r); 
            lo = P.Delta_Low(r); hi = P.Delta_High(r);
            eff = sprintf('Δ = %s ms [%s, %s]', sgn(d), fnum(lo, 2), fnum(hi, 2));
        case 'coverage'
            rr = P.OR(r); 
            lo = P.OR_Low(r); hi = P.OR_High(r);
            eff = sprintf('OR = %s [%s, %s]', fnum(rr,2), fnum(lo, 2), fnum(hi, 2));
        otherwise
            rr = P.RR(r); 
            lo = P.RR_Low(r); hi = P.RR_High(r);
            eff = sprintf('RR = %s [%s, %s]', fnum(rr,2), fnum(lo, 2), fnum(hi, 2));
    end

    dirTxt = direction_text(E1, E2, L1, L2, u);
    fprintf('  %s vs %s: %s; EMMs: %s %s (%s) vs %s %s (%s); p = %s (pFDR = %s)\n', L1, L2, eff, fnum(E1,2), u, L1, fnum(E2,2), u, L2, fp(P.pValue(r)), fp(P.pFDR(r)));

    if ~isempty(dirTxt)
        fprintf('    Direction: %s\n', dirTxt);
    end
    disp(' ');
end
end

function print_pairs_simple(P, T, measure, withinF, byF, alphaFDR)

if ~ismember('pFDR', P.Properties.VariableNames)
    P.pFDR = bh_adjust(P.pValue);
end

sig = P.pFDR < alphaFDR & ~isnan(P.pFDR);

if ~any(sig)
    fprintf('  No pairwise survived FDR.\n');
    return;
end

u = units_for(measure);
rows = find(sig);

for r = rows(:).'
    ctx = string(P.(byF)(r));
    L1  = string(P.Level1(r)); 
    L2 = string(P.Level2(r));

    % EMMs from table that contains both factors

    if strcmpi(withinF,'Group')
        E1 = pick_emm2(T, 'Group', L1, 'Condition', ctx);
        E2 = pick_emm2(T, 'Group', L2, 'Condition', ctx);
    else
        E1 = pick_emm2(T, 'Condition', L1, 'Group', ctx);
        E2 = pick_emm2(T, 'Condition', L2, 'Group', ctx);
    end

    switch lower(measure)
        case 'duration'
            d  = P.Delta(r);
            lo = P.Delta_Low(r);
            hi = P.Delta_High(r);
            eff = sprintf('Δ = %s ms [%s, %s]', sgn(d), fnum(lo,2), fnum(hi,2));
        case 'coverage'
            or = P.OR(r);
            lo = P.OR_Low(r); 
            hi = P.OR_High(r);
            eff = sprintf('OR = %s [%s, %s]', fnum(or,2), fnum(lo,2), fnum(hi,2));
        otherwise
            rr = P.RR(r);
            lo = P.RR_Low(r);
            hi = P.RR_High(r);
            eff = sprintf('RR = %s [%s, %s]', fnum(rr,2), fnum(lo,2), fnum(hi,2));
    end

    dirTxt = direction_text(E1, E2, L1, L2, u);

    fprintf('  %s vs %s @ %s=%s: %s; EMMs: %s %s (%s) vs %s %s (%s); p = %s (pFDR = %s)\n', ...
        L1, L2, byF, ctx, eff, fnum(E1, 2), u, L1, fnum(E2, 2), u, L2, fp(P.pValue(r)), fp(P.pFDR(r)));

    if ~isempty(dirTxt)
        fprintf('    Direction: %s\n', dirTxt);
    end
end
end

function x = pick_emm(T, fac, lvl)

mask = T.(fac) == categorical(lvl, categories(T.(fac)));
x = NaN;

if any(mask)
    x = T.EMM_BT(find(mask, 1, 'first')); 
end
end

function x = pick_emm2(T, f1, l1, f2, l2)

mask = T.(f1) == categorical(l1, categories(T.(f1))) & T.(f2) == categorical(l2, categories(T.(f2)));
x = NaN;

if any(mask)
    x = T.EMM_BT(find(mask, 1, 'first')); 
end
end

function val = get_nested(S, path)

val = NaN;

try
    for i = 1: numel(path)
        S = S.(path{i});
    end
    val = S;
catch
end
end

function u = units_for(measure)

switch measure
    case 'Duration'
        u = 'ms';
    case 'Coverage'
        u = '%';
    otherwise
        u = 'Hz';
end
end

function s = direction_text(x1, x2, l1, l2, u)

if isnan(x1) || isnan(x2)
    s = ''; 
    return; 
end

if x1 > x2
    s = sprintf('%s > %s by %s %s', l1, l2, fnum(x1-x2,2), u);
elseif x2 > x1
    s = sprintf('%s > %s by %s %s', l2, l1, fnum(x2-x1,2), u);
else
    s = sprintf('%s = %s', l1, l2);
end
end

function s = fp(p)

if isnan(p)
    s = 'NA'; 
elseif p < 1e-4
    s = '<1e-4';
else
    s = sprintf('%.3f', p);
end
end

function s = fnum(x, d)

if isnan(x)
    s = 'NA';
else
    s = sprintf(['%0.', num2str(d), 'f'], x);
end
end

function s = fdf(x)

if isnan(x)
    s = 'NA'; 
else
    s = sprintf('%.1f', x);
end
end

function s = sgn(x)

if isnan(x)

    s = 'NA';

elseif x >= 0
    s = sprintf('+%.2f', x);
else
    s = sprintf('%.2f', x);
end
end

function s = star(tf)

if tf 
    s = '★';
else
    s = ''; 
end
end

function p_adj = bh_adjust(p)

p = p(:); m = numel(p);
[ps, idx] = sort(p);
q = ps .* m ./ (1:m)';

for i = m - 1: -1: 1
    q(i) = min(q(i), q(i + 1));
end

p_adj = nan(size(p)); p_adj(idx) = q;

end
