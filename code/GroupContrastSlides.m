function GroupContrastSlides(CompStats, erp, Conds, Groups, varargin) % Make one figure per Target group (≠ RefGroup), with 3 circular plots (one per Condition) % for a chosen ERP. Recomputes BH–FDR within each (ERP × Condition × TargetGroup) across 42 pairs.

p = inputParser;

addParameter(p, 'RefGroup', ''); % auto-detect if empty
addParameter(p, 'GroupsToPlot', {}); % auto = all except RefGroup
addParameter(p, 'Alpha', []);
addParameter(p, 'MicrostateLabels', []);
addParameter(p, 'NodeRadius', []);
addParameter(p, 'ArrowOffset', []);
addParameter(p, 'WidthDomain', []);
addParameter(p, 'CurveMagnitude', []);
addParameter(p, 'MinWidth', []);
addParameter(p, 'MaxWidth', []);
addParameter(p, 'GreyColor', []);
 addParameter(p, 'ShaftGap', []);
 addParameter(p, 'AboveColor', []);
 addParameter(p, 'BelowColor', []);
 addParameter(p, 'GreyLineWidth', []); 

parse(p, varargin{:}); 
opt = p.Results;

% detect reference group

refGroup = char(opt.RefGroup);

if isempty(refGroup)
    refGroup = detectRefGroup(CompStats, erp, Conds, Groups);
end

if isempty(refGroup)
    error('Could not detect RefGroup from CompStats; please pass ''RefGroup'' explicitly.');
end

% groups to plot

if isempty(opt.GroupsToPlot)
    groups_to_plot = Groups(~strcmpi(Groups, refGroup));
else
    groups_to_plot = cellstr(opt.GroupsToPlot);
    groups_to_plot = groups_to_plot(~strcmpi(groups_to_plot, refGroup) & ismember(groups_to_plot, Groups));
end

if isempty(groups_to_plot)
    warning('No target groups to plot (all are RefGroup=%s). Nothing to do.', refGroup);
    return;
end

nCols = numel(Conds);

for gi = 1: numel(groups_to_plot)

    grp = groups_to_plot{gi};
    figName = sprintf('%s vs %s — %s', strrep(grp, '_', ' '), strrep(refGroup, '_', ' '), erp);
    figure('Color', 'w', 'Units', 'normalized', 'Position', [0.0500 0.6022 0.2750 0.1489], 'Name', figName, 'Visible', 'on');
    tl = tiledlayout(1, nCols, 'Padding', 'none', 'TileSpacing', 'none');

    for ci = 1: nCols
        cond = Conds{ci};
        nexttile(tl, ci);

        % recompute per-group FDR (preferred) or fall back to stored edges

        E = table([], [], [], [], [], [], 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'qFDR'});

        if isfield(CompStats, erp) && isfield(CompStats.(erp), cond) && isfield(CompStats.(erp).(cond), grp)
            node = CompStats.(erp).(cond).(grp);
            can_recompute = isfield(node, 'p') && ~isempty(node.p) && isfield(node, 'IRR') && ~isempty(node.IRR);

            if can_recompute
                pMat = node.p;
                IRR = node.IRR;
                have_CI = isfield(node, 'IRR_CIlo') && isfield(node, 'IRR_CIhi');

                if have_CI
                    CIlo = node.IRR_CIlo;
                    CIhi = node.IRR_CIhi;
                else
                    CIlo = NaN(7);
                    CIhi = NaN(7);
                end

                offdiag = ~eye(7);
                p_lin = pMat(offdiag);
                keep = ~isnan(p_lin);
                qMat = NaN(7);

                if any(keep)
                    q_all = bh_adjust(p_lin(keep));
                    q_only = NaN(size(p_lin));
                    q_only(keep) = q_all;
                    qMat(offdiag) = q_only;
                end

                sig = qMat < opt.Alpha;
                [sSig, tSig] = find(sig);

                if ~isempty(sSig)
                    irr_vec = IRR(sub2ind([7, 7], sSig, tSig));
                    lo_vec = CIlo(sub2ind([7, 7], sSig, tSig));
                    hi_vec = CIhi(sub2ind([7, 7], sSig, tSig));
                    q_vec = qMat(sub2ind([7, 7], sSig, tSig));
                    E = table(sSig, tSig, irr_vec, lo_vec, hi_vec, q_vec, 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'qFDR'});
                    E.srcLabel = reshape(string(opt.MicrostateLabels(E.src)), [], 1);
                    E.tgtLabel = reshape(string(opt.MicrostateLabels(E.tgt)), [], 1);
                    dlab = strings(height(E), 1);
                    dlab(E.IRR > 1) = "Higher vs " + string(refGroup); dlab(E.IRR < 1) = "Lower vs " + string(refGroup);
                    E.Direction = dlab;
                else
                    E = table([], [], [], [], [], [], strings(0, 1), 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'qFDR', 'Direction'});
                end

            elseif isfield(node, 'edges')
                E = node.edges;
                if ~ismember('Direction', E.Properties.VariableNames)
                    if ~isempty(E)
                        E.Direction = strings(height(E), 1);
                        E.Direction(E.IRR > 1) = "Higher vs " + string(refGroup); E.Direction(E.IRR < 1) = "Lower vs " + string(refGroup);
                    else
                        E = table([], [], [], [], [], [], strings(0, 1), 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'qFDR', 'Direction'});
                    end
                end
            end
        end

        plotCircularMicrostateTransitions(E, 'Labels', opt.MicrostateLabels, 'NodeRadius', opt.NodeRadius, 'ArrowOffset', opt.ArrowOffset, 'WidthDomain', opt.WidthDomain, 'MinWidth', opt.MinWidth, 'MaxWidth', opt.MaxWidth, 'CurveMagnitude', opt.CurveMagnitude, 'GreyMesh', true, 'GreyLineWidth', opt.GreyLineWidth, 'AboveColor', opt.AboveColor, 'BelowColor', opt.BelowColor, 'GreyColor', opt.GreyColor, 'ShaftGap', opt.ShaftGap);
    end

    drawnow;
end
end

function ref = detectRefGroup(CompStats, erp, Conds, Groups)
ref = '';

if ~isfield(CompStats, erp)
    return;
end

S1 = CompStats.(erp);

for c = 1: numel(Conds)
    cn = Conds{c};

    if ~isfield(S1, cn)
        continue;
    end

    S2 = S1.(cn);
    gnames = fieldnames(S2);

    for k = 1: numel(gnames)
        g = gnames{k};
        if any(strcmpi(g, Groups)) && isfield(S2.(g), 'RefGroup')
            ref = char(S2.(g).RefGroup);
            if ~isempty(ref)
                return;
            end
        end
    end
end
end

function q = bh_adjust(p)
p = p(:); [ps, ord] = sort(p, 'ascend');
m = numel(ps);
ranks = (1: m)';
qtmp = (m ./ ranks) .* ps;
qtmp = flipud(cummin(flipud(qtmp)));
q = NaN(size(ps));
q(ord) = min(1, qtmp);
end
