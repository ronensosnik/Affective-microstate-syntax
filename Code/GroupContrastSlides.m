function GroupContrastSlides(CompStats, erp, Conds, Groups, varargin) % Make one figure per Target group (≠ RefGroup), with 3 circular plots (one per Condition) % for a chosen ERP. Recomputes BH–FDR within each (ERP × Condition × TargetGroup) across 42 pairs.

p = inputParser;

% ====================================================================%

% Node radius and label

addParameter(p, 'NodeRadius', []);
addParameter(p, 'NodeLabels', []);

% Arrow width, color, starting and ending points

addParameter(p, 'WidthDomain', []);
addParameter(p, 'MinArrowWidth', []);
addParameter(p, 'MaxArrowWidth', []);

addParameter(p, 'GreyColor', []);
addParameter(p, 'GreyLineWidth', []);
addParameter(p, 'AboveColor', []);
addParameter(p, 'BelowColor', []);

% Bubbles radii, color, and edge

addParameter(p, 'DurRadiusMin', []);
addParameter(p, 'DurRadiusMax', []);
addParameter(p, 'DurDomain', []);
addParameter(p, 'DurFaceColor', []);
addParameter(p, 'DurEdgeColor', []);
addParameter(p, 'DurEdgeWidth', []);

% ====================================================================%

parse(p, varargin{:});
opt = p.Results;

% detect reference group

refGroup = 'HC';
groups_to_plot = Groups(~strcmpi(Groups, 'HC'));

nCols = numel(Conds);

for gi = 1: numel(groups_to_plot)

    grp = groups_to_plot{gi};
    figName = sprintf('%s vs %s — %s', strrep(grp, '_', ' '), strrep('HC', '_', ' '), erp);
    figure('Color', 'w', 'Units', 'normalized', 'Position', [0.0500 0.6022 0.2750 0.1489], 'Name', figName, 'Visible', 'on');
    tl = tiledlayout(1, nCols, 'Padding', 'none', 'TileSpacing', 'none');

    for ci = 1: nCols

        cond = Conds{ci};
        nexttile(tl, ci);

        % recompute per-group FDR

        E = table([], [], [], [], [], [], 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'pFDR'});

        if isfield(CompStats, erp) && isfield(CompStats.(erp), cond) && isfield(CompStats.(erp).(cond), grp)
            node = CompStats.(erp).(cond).(grp);
            can_recompute = isfield(node, 'p') && ~isempty(node.p) && isfield(node, 'IRR') && ~isempty(node.IRR);

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
            pFDR_Mat = NaN(7);

            if any(keep)
                pFDR_all = bh_adjust(p_lin(keep));
                pFDR_only = NaN(size(p_lin));
                pFDR_only(keep) = pFDR_all;
                pFDR_Mat(offdiag) = pFDR_only;
            end

            sig = pFDR_Mat < 0.05;
            [sSig, tSig] = find(sig);

            if ~isempty(sSig)
                irr_vec = IRR(sub2ind([7, 7], sSig, tSig));
                lo_vec = CIlo(sub2ind([7, 7], sSig, tSig));
                hi_vec = CIhi(sub2ind([7, 7], sSig, tSig));
                pFDR_vec = pFDR_Mat(sub2ind([7, 7], sSig, tSig));
                E = table(sSig, tSig, irr_vec, lo_vec, hi_vec, pFDR_vec, 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'pFDR'});
                E.srcLabel = reshape(string(opt.NodeLabels(E.src)), [], 1);
                E.tgtLabel = reshape(string(opt.NodeLabels(E.tgt)), [], 1);
                dlab = strings(height(E), 1);
                dlab(E.IRR > 1) = "Higher vs HC";
                dlab(E.IRR < 1) = "Lower vs HC";
                E.Direction = dlab;
                E.source_dur = node.edges.source_dur;
                E.target_dur = node.edges.target_dur;
            else
                E = table([], [], [], [], [], [], strings(0, 1), [], [], 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'ObsOverExp', 'Direction', 'source_dur', 'target_dur'});
            end

        end

        plotCircularMicrostateTransitions(E, 'NodeLabels', opt.NodeLabels, 'NodeRadius', opt.NodeRadius, 'GreyMesh', 'true', 'GreyLineWidth', opt.GreyLineWidth, 'GreyColor', opt.GreyColor, 'WidthDomain', opt.WidthDomain, 'MinArrowWidth', opt.MinArrowWidth, 'MaxArrowWidth', opt.MaxArrowWidth, 'AboveColor', opt.AboveColor, 'BelowColor', opt.BelowColor, 'DurDomain', opt.DurDomain, 'DurRadiusMin', opt.DurRadiusMin, 'DurRadiusMax', opt.DurRadiusMax, 'DurFaceColor', opt.DurFaceColor, 'DurEdgeColor', opt.DurEdgeColor, 'DurEdgeWidth', opt.DurEdgeWidth);
    end

    drawnow;
end
end

function pFDR = bh_adjust(p)
p = p(:);
[ps, ord] = sort(p, 'ascend');
m = numel(ps);
ranks = (1: m)';
qtmp = (m ./ ranks) .* ps;
qtmp = flipud(cummin(flipud(qtmp)));
pFDR = NaN(size(ps));
pFDR(ord) = min(1, qtmp);
end