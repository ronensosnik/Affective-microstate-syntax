function [overall, byERP] = summarize_q1_robustness(Base, Med, ERPs, Groups, Conds)

% Compare Med to Base on baseline significant edges: direction retention and significance retention.

total_base = 0;
retain_dir = 0;
retain_sig = 0;
per_erp = struct();

for e = 1: numel(ERPs)
    erp = ERPs{e};
    bE = 0;
    dE = 0;
    sE = 0;

    for g = 1:numel(Groups)
        grp = Groups{g};

        for c = 1:numel(Conds)
            cond = Conds{c};

            if ~isfield(Base.(erp).(grp), cond) || ~isfield(Med.(erp).(grp), cond)
                continue;
            end

            B = Base.(erp).(grp).(cond); M = Med.(erp).(grp).(cond);
            sigB = B.sig;
            if isempty(sigB)
                continue;
            end

            for s = 1:7
                for t = 1:7

                    if s==t
                        continue;
                    end

                    if sigB(s,t)
                        total_base = total_base + 1;
                        bE = bE + 1;
                        dirB = sign(B.IRR(s,t) - 1);
                        dirM = sign(M.IRR(s,t) - 1);
                        if ~isnan(dirB) && ~isnan(dirM) && dirB==dirM
                            retain_dir = retain_dir + 1; dE = dE + 1;
                        end
                        if isfield(M,'sig') && ~isempty(M.sig) && M.sig(s,t)
                            retain_sig = retain_sig + 1; sE = sE + 1;
                        end
                    end
                end
            end
        end
    end
    per_erp.(erp) = struct('BaselineEdges', bE, 'DirectionRetained', dE, 'SignificanceRetained', sE, 'DirPct', pct(dE,bE), 'SigPct', pct(sE,bE));
end

overall = struct('BaselineEdges', total_base, 'DirectionRetained', retain_dir, 'SignificanceRetained', retain_sig, 'DirPct', pct(retain_dir,total_base), 'SigPct', pct(retain_sig,total_base));

erp_rows = table('Size', [numel(ERPs) 5], 'VariableTypes', {'string', 'double', 'double', 'double', 'string'}, 'VariableNames', {'ERP', 'BaselineEdges', 'DirectionRetained', 'SignificanceRetained', 'Percentages'});

for e = 1: numel(ERPs)
    S = per_erp.(ERPs{e});
    erp_rows.ERP(e) = string(ERPs{e});
    erp_rows.BaselineEdges(e) = S.BaselineEdges;
    erp_rows.DirectionRetained(e) = S.DirectionRetained;
    erp_rows.SignificanceRetained(e) = S.SignificanceRetained;
    erp_rows.Percentages(e) = sprintf('Dir=%s, Sig=%s', pct(S.DirectionRetained,S.BaselineEdges), pct(S.SignificanceRetained,S.BaselineEdges));
end
byERP = erp_rows;
end