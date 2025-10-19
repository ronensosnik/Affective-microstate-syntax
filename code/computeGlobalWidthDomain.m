function WidthDomain = computeGlobalWidthDomain(TransStats_cellFDR, CompStats, ERPs, Conds, Groups)

% Returns [0, cap], where cap = max |log(IRR)| across:
% Q1: TransStats_cellFDR.(ERP).(Group).(Cond).IRR
% Q2: CompStats.(ERP).(Cond).(Target).IRR
% Off-diagonal pairs only; ignores NaNs/infs.

vals = [];

%%% Q1 %%%%

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

            if ~isfield(node, 'IRR') || isempty(node.IRR)
                continue;
            end

            IRR = node.IRR;
            mask = ~eye(7) & isfinite(IRR);
            v = abs(log(IRR(mask)));
            vals = [vals; v(:)];
        end
    end
end

%%% Q2 %%%%

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

            if ~isfield(node, 'IRR') || isempty(node.IRR) || ~ismember(tgt, Groups)
                continue;
            end

            IRR = node.IRR;
            mask = ~eye(7) & isfinite(IRR);
            v = abs(log(IRR(mask)));
            vals = [vals; v(:)];
        end
    end
end

vals = vals(isfinite(vals) & vals >= 0);

if isempty(vals)
    cap = 1;
else
    cap = max(vals);
    if cap <= 0
        cap = 1;
    end
end

WidthDomain = [0 cap];
end
