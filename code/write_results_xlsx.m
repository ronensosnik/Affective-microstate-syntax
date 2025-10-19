function write_results_xlsx(save_path, ERPs, Conds, Groups, TransStats_cellFDR, CompStats, RefGroup)

% Write only significant Q1/Q2 results (q < 0.05) to Excel, one file per ERP.
%
% Sheets:
%   - Q1_independence
%   - Q2_vs_<RefGroup>

    ALPHA_Q = 0.05;

    for e = 1: numel(ERPs)
        erp = ERPs{e};

        % ------- Q1 table (within-cell independence models) -------

        Tq1 = buildQ1TableForERP(erp, Conds, Groups, TransStats_cellFDR, ALPHA_Q);

        % ------- Q2 table (anchored contrasts vs reference) -------

        Tq2 = buildQ2TableForERP(erp, Conds, Groups, CompStats, RefGroup, ALPHA_Q);

        % ------- write Excel (one file per ERP, two sheets) -------

        xname = fullfile(save_path, sprintf('Microstates_%s_Q1Q2.xlsx', erp));

        if exist(xname, 'file')
            delete(xname);
        end

        if ~isempty(Tq1)
            writetable(Tq1, xname, 'Sheet', 'Q1_independence');
        else
            writetable(cell2table({"No Q1 rows (q < 0.05) for this ERP"}, 'VariableNames', {'Note'}), xname, 'Sheet', 'Q1_independence');
        end

        sheetQ2 = sprintf('Q2_vs_%s', strrep(RefGroup, '_', ' '));

        if ~isempty(Tq2)
            writetable(Tq2, xname, 'Sheet', sheetQ2);
        else
            writetable(cell2table({"No Q2 rows (q < 0.05) for this ERP"}, 'VariableNames', {'Note'}), xname, 'Sheet', sheetQ2);
        end
    end
end

function T = buildQ1TableForERP(erp, Conds, Groups, TransStats, alpha_q)

% Build Q1 table for a single ERP and keep only rows with q < alpha_q.
% Also convert src/tgt 1..7 to letters A..G.

    rows = {};
    Labels = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};

    for gi = 1: numel(Groups)
        grp = Groups{gi};

        if ~isfield(TransStats, erp) || ~isfield(TransStats.(erp), grp)
            continue;
        end

        for ci = 1: numel(Conds)
            cond = Conds{ci};

            if ~isfield(TransStats.(erp).(grp), cond)
                continue;
            end

            node = TransStats.(erp).(grp).(cond);

            if ~isfield(node, 'IRR')
                node.IRR = NaN(7);
            end

            if ~isfield(node, 'p')
                node.p = NaN(7);
            end

            if ~isfield(node, 'q')
                node.q = NaN(7);
            end

            if ~isfield(node, 'nRows')
                node.nRows = zeros(7);
            end

            if ~isfield(node, 'ObsOverExp')
                node.ObsOverExp = NaN(7);
            end

            % Reconstruct Wald 95% CIs for IRR if not stored

            if ~isfield(node, 'IRR_CIlo') || ~isfield(node, 'IRR_CIhi') || isempty(node.IRR_CIlo) || isempty(node.IRR_CIhi)
                if isfield(node, 'beta0') && isfield(node, 'se') && ~isempty(node.beta0) && ~isempty(node.se)
                    node.IRR_CIlo = exp(node.beta0 - 1.96 .* node.se);
                    node.IRR_CIhi = exp(node.beta0 + 1.96 .* node.se);
                else
                    node.IRR_CIlo = NaN(7);
                    node.IRR_CIhi = NaN(7);
                end
            end

            offdiag = find(~eye(7));
            [S, Tt] = ind2sub([7, 7], offdiag);

            for k = 1: numel(offdiag)
                s = S(k);
                t = Tt(k);

                rows(end + 1, :) = {erp, grp, cond, s, t, safeGet(node.IRR, s, t), safeGet(node.IRR_CIlo, s, t), safeGet(node.IRR_CIhi, s, t), safeGet(node.p, s, t), safeGet(node.q, s, t), safeGet(node.nRows, s, t), safeGet(node.ObsOverExp, s, t)                 };
            end
        end
    end

    if isempty(rows)
        T = table();
    else
        T = cell2table(rows, 'VariableNames', {'ERP', 'Group','Condition', 'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'p', 'q', 'nRows', 'ObsOverExp'});

        % --- keep only significant rows

        mask = isfinite(T.q) & T.q < alpha_q;
        T = T(mask, :);

        % --- convert src/tgt numbers to letters

        if ~isempty(T)
            idx = T.src;
            src_labels = string(Labels(idx));

            if isrow(src_labels)
                src_labels = src_labels.';
            end

            T.src = src_labels;

            idx = T.tgt;
            tgt_labels = string(Labels(idx));

            if isrow(tgt_labels)
                tgt_labels = tgt_labels.';
            end
            
            T.tgt = tgt_labels;

        end
    end
end

function T = buildQ2TableForERP(erp, Conds, Groups, CompStats, RefGroup, alpha_q)

% Build Q2 table for a single ERP vs RefGroup and keep only rows with q < alpha_q.
% Also convert src/tgt 1..7 to letters A..G.

    rows = {};
    Labels = {'A', 'B', 'C', 'D', 'E', 'F', 'G'};

    if ~isfield(CompStats, erp)
        T = table();
        return;
    end

    for ci = 1: numel(Conds)
        cond = Conds{ci};

        if ~isfield(CompStats.(erp), cond)
            continue;
        end

        S = CompStats.(erp).(cond);

        tgNames = fieldnames(S);
        tgNames = tgNames(ismember(tgNames, Groups) & ~strcmpi(tgNames, RefGroup));

        for ti = 1: numel(tgNames)

            tgt = tgNames{ti};
            node = S.(tgt);

            if ~isfield(node, 'IRR')
                node.IRR = NaN(7);
            end

            if ~isfield(node, 'IRR_CIlo')
                node.IRR_CIlo = NaN(7);
            end

            if ~isfield(node, 'IRR_CIhi')
                node.IRR_CIhi = NaN(7);
            end

            if ~isfield(node, 'p')
                node.p = NaN(7);
            end

            if ~isfield(node, 'q')
                node.q = NaN(7);
            end

            offdiag = find(~eye(7));
            [Sidx, Tidx] = ind2sub([7, 7], offdiag);

            for k = 1: numel(offdiag)

                s = Sidx(k);
                t = Tidx(k);

                rows(end + 1, :) = {erp, cond, RefGroup, tgt, s, t, safeGet(node.IRR, s, t), safeGet(node.IRR_CIlo, s, t), safeGet(node.IRR_CIhi, s, t), safeGet(node.p, s, t), safeGet(node.q, s, t)}; 
            end
        end
    end

    if isempty(rows)
        T = table();
    else

        T = cell2table(rows, 'VariableNames', {'ERP', 'Condition', 'RefGroup', 'TargetGroup', 'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'p', 'q'});

        % --- keep only significant rows

        mask = isfinite(T.q) & T.q < alpha_q;
        T = T(mask, :);

        % --- convert src/tgt numbers to letters

        if ~isempty(T)
            idx = T.src;
            src_labels = string(Labels(idx));

            if isrow(src_labels)
                src_labels = src_labels.';
            end

            T.src = src_labels;

            idx = T.tgt;
            tgt_labels = string(Labels(idx));

            if isrow(tgt_labels)
                tgt_labels = tgt_labels.';
            end

            T.tgt = tgt_labels;
        end
    end
end

function v = safeGet(M, s, t)
    try
        v = M(s, t);
    catch
        v = NaN;
    end
end