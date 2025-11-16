function [overall, byERP] = summarize_q2_robustness(Base, Med, ERPs, Conds, Groups)

% Compare Med vs Base on baseline significant edges per Target×ERP×Condition.
% Counts: total baseline edges; of those, retain direction? retain significance under Med?

total_base = 0;
retain_dir = 0;
retain_sig = 0;
per_erp = struct();

for e = 1: numel(ERPs)
    erp = ERPs{e};
    bE = 0;
    dE = 0;
    sE = 0;

    for c = 1: numel(Conds)

        cond = Conds{c};

        if ~isfield(Base.(erp),cond) || ~isfield(Med.(erp),cond)
            continue;
        end

        tgs = fieldnames(Base.(erp).(cond));

        % keep only true target groups

        tgs = tgs(ismember(tgs, Groups) & ~strcmp(tgs, 'FDR_SCOPE'));

        for k = 1: numel(tgs)
            tgt = tgs{k};

            B = Base.(erp).(cond).(tgt);

            if ~isfield(Med.(erp).(cond), tgt)
                continue;
            end

            M = Med.(erp).(cond).(tgt);

            sigB = B.sig;

            if isempty(sigB)
                continue;
            end

            for s = 1: 7
                for t = 1: 7

                    if s == t
                        continue;
                    end

                    if sigB(s, t)
                        total_base = total_base + 1;
                        bE = bE + 1;
                        dirB = sign(B.IRR(s, t) - 1);
                        dirM = sign(M.IRR(s, t) - 1);

                        if ~isnan(dirB) && ~isnan(dirM) && dirB == dirM
                            retain_dir = retain_dir + 1;
                            dE = dE + 1;
                        end

                        if isfield(M, 'sig') && ~isempty(M.sig) && M.sig(s, t)
                            retain_sig = retain_sig + 1;
                            sE = sE + 1;
                        end
                    end
                end
            end
        end
    end
    per_erp.(erp) = struct('BaselineEdges', bE, 'DirectionRetained', dE, 'SignificanceRetained', sE, 'DirPct', pct(dE,bE), 'SigPct', pct(sE, bE));
end

overall = struct('BaselineEdges', total_base, 'DirectionRetained', retain_dir, 'SignificanceRetained', retain_sig, 'DirPct', pct(retain_dir, total_base), 'SigPct', pct(retain_sig,total_base));

% Make a by-ERP table

byERP = table('Size',[numel(ERPs) 5], 'VariableTypes', {'string', 'double', 'double', 'double', 'string'}, 'VariableNames', {'ERP', 'BaselineEdges', 'DirectionRetained', 'SignificanceRetained', 'Percentages'});

for e = 1: numel(ERPs)
    S = per_erp.(ERPs{e});
    byERP.ERP(e) = string(ERPs{e});
    byERP.BaselineEdges(e) = S.BaselineEdges;
    byERP.DirectionRetained(e) = S.DirectionRetained;
    byERP.SignificanceRetained(e) = S.SignificanceRetained;
    byERP.Percentages(e) = sprintf('Dir=%s, Sig=%s', pct(S.DirectionRetained, S.BaselineEdges), pct(S.SignificanceRetained, S.BaselineEdges));
end

end