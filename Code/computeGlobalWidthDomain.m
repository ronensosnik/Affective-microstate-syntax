function [WidthDomain, DurDomain] = computeGlobalWidthDomain(TransStats_cellFDR, CompStats, ERPs, Conds, Groups)

% WidthDomain: [cap(1) cap(2)], where cap = [min |log(IRR)| max |log(IRR)|] across Q1 + Q2 (off-diagonal).
% DurDomain: [dmin dmax], where dmin/dmax are taken from duration columns on significant edges across Q1 + Q2.

vals = [];
durVals = [];

% --------------------------
% Q1 width + duration
% --------------------------

for e = 1: numel(ERPs)

    erp = ERPs{e};

    if ~isfield(TransStats_cellFDR, erp)
        continue;
    end

    for g = 1: numel(Groups)

        grp = Groups{g};

        if ~isfield(TransStats_cellFDR.(erp), grp)
            continue;
        end

        for c = 1: numel(Conds)

            cond = Conds{c};

            if ~isfield(TransStats_cellFDR.(erp).(grp), cond)
                continue;
            end

            node = TransStats_cellFDR.(erp).(grp).(cond);

            if isfield(node, 'edges') && ~isempty(node.edges)
                clear IRR v;
                IRR = node.edges.IRR;
                v = abs(log(IRR));
                vals = [vals; v(:)]; % transition strength
                durVals = [durVals; collectDurationsFromEdges(node.edges)]; % duration spent at source and target
            end
        end
    end
end

% --------------------------
% Q2 width + duration
% --------------------------

for e = 1: numel(ERPs)

    erp = ERPs{e};

    if ~isfield(CompStats, erp)
        continue;
    end

    for c = 1: numel(Conds)

        cond = Conds{c};

        if ~isfield(CompStats.(erp), cond)
            continue;
        end

        tgNames = fieldnames(CompStats.(erp).(cond));

        for k = 1: numel(tgNames)

            tgt = tgNames{k};
            node = CompStats.(erp).(cond).(tgt);

            if ~ismember(tgt, Groups)
                continue;
            end

            if isfield(node, 'edges') && ~isempty(node.edges)
                clear IRR v;
                IRR = node.edges.IRR;
                v = abs(log(IRR));
                vals = [vals; v(:)]; % transition strength
                durVals = [durVals; collectDurationsFromEdges(node.edges)]; % duration spent at source and target
            end
        end
    end
end

vals = vals(isfinite(vals) & vals >= 0);

WidthDomain = [min(vals) max(vals)];

durVals = durVals(isfinite(durVals));

dmin = min(durVals);
dmax = max(durVals);

if ~isfinite(dmin) || ~isfinite(dmax) || dmax <= dmin
    DurDomain = [dmin max(dmin + 1, dmax)];
else
    DurDomain = [dmin dmax];
end
end

function dv = collectDurationsFromEdges(E)

dv = [];

if ~istable(E)
    return;
end

cols = E.Properties.VariableNames;

srcCol = pickFirstExisting(cols, {'source_dur', 'Dur_src', 'Dur_source'});
tgtCol = pickFirstExisting(cols, {'target_dur', 'Dur_trg', 'Dur_target'});

if ~isempty(srcCol)
    x = E.(srcCol);
    dv = [dv; x(:)];
end

if ~isempty(tgtCol)
    x = E.(tgtCol);
    dv = [dv; x(:)];
end

end

function name = pickFirstExisting(cols, candidates)

name = '';

for i = 1: numel(candidates)
    if any(strcmp(cols, candidates{i}))
        name = candidates{i};
        return;
    end
end

end