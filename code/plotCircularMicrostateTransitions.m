function plotCircularMicrostateTransitions(edges, varargin)

% Circular microstate graph for one condition (curved, filled arrowheads) with duration bubbles for significant edges.
%
% Options:
% 'NodeLabels' -> node labels (7)
% 'NodeRadius' -> circle radius for nodes
% 'WidthDomain' -> [cap(1) cap(2)] maps |log(IRR)| globally to [MinArrowWidth, MaxArrowWidth]
% 'MinArrowWidth', 'MaxArrowWidth'
% 'GreyMesh', 'GreyLineWidth', 'GreyColor'
% 'AboveColor', 'BelowColor'

% Duration bubbles (for significant edges only):
% Reads durations from edges table columns:
% source_dur / target_dur (preferred)
% 'DurDomain' -> [min max] for duration scaling across all Q1/Q2 plots (recommended)
% 'DurRadiusMin', 'DurRadiusMax'
% 'DurFaceColor', 'DurEdgeColor', 'DurEdgeWidth'

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

% Mesh

addParameter(p, 'GreyMesh', []);

% ====================================================================%

parse(p, varargin{:});
opt = p.Results;

labels = opt.NodeLabels;

if isstring(labels)
    labels = cellstr(labels);
end

% Normalize / fill edge table

if isempty(edges)
    edges = table([], [], [], [], [], [], strings(0, 1), [], [], 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'ObsOverExp', 'Direction', 'source_dur', 'target_dur'});
end

if ismember('Source', edges.Properties.VariableNames)
    edges.Properties.VariableNames{'Source'} = 'src';
end

if ismember('Target', edges.Properties.VariableNames)
    edges.Properties.VariableNames{'Target'} = 'tgt';
end

if ~ismember('Direction', edges.Properties.VariableNames)
    edges.Direction = strings(height(edges), 1);
    edges.Direction(edges.IRR > 1) = "Above expected";
    edges.Direction(edges.IRR < 1) = "Below expected";
end

cla;
hold on;
axis equal off;

pad = 1.35;
xlim([-pad pad]);
ylim([-pad pad]);

% Node positions: A at 0° (right), others CCW

theta = 2 * pi * (0: 6) / 7;
x = cos(theta);
y = sin(theta);

nodeR = opt.NodeRadius;

% Background mesh

if opt.GreyMesh
    for i = 1: 7
        for j = i + 1: 7
            plot([x(i) x(j)], [y(i) y(j)], ':', 'Color', opt.GreyColor, 'LineWidth', opt.GreyLineWidth);
        end
    end
end

% Collect significant edges into matrices

sig = false(7);
irr = NaN(7);
dirSign = zeros(7);

for r = 1: height(edges)
    i = edges.src(r);
    j = edges.tgt(r);

    if i == j || any([i j] < 1) || any([i j] > 7)
        continue;
    end

    if ~isfinite(edges.IRR(r)) || edges.IRR(r) <= 0
        continue;
    end

    sig(i, j) = true;
    irr(i, j) = edges.IRR(r);

    d = lower(string(edges.Direction(r)));

    if startsWith(d, "above") || startsWith(d, "higher")
        dirSign(i, j) = +1;
    elseif startsWith(d, "below") || startsWith(d, "lower")
        dirSign(i, j) = -1;
    else
        dirSign(i, j) = sign(log(max(eps, edges.IRR(r))));
    end
end

% Width mapping

useDomain = ~isempty(opt.WidthDomain) && numel(opt.WidthDomain) == 2 && diff(opt.WidthDomain) > 0;
wmin = opt.MinArrowWidth;
wmax = opt.MaxArrowWidth;

if useDomain
    cap = opt.WidthDomain;
else
    cap = [];
end

mapEffectToWidth = @(val) local(val, cap, wmin, wmax);

% Duration columns (support both naming conventions)

srcDurCol = '';
tgtDurCol = '';

if ismember('source_dur', edges.Properties.VariableNames)
    srcDurCol = 'source_dur';
elseif ismember('Dur_src', edges.Properties.VariableNames)
    srcDurCol = 'Dur_src';
end

if ismember('target_dur', edges.Properties.VariableNames)
    tgtDurCol = 'target_dur';
elseif ismember('Dur_trg', edges.Properties.VariableNames)
    tgtDurCol = 'Dur_trg';
end

hasDur = ~isempty(srcDurCol) && ~isempty(tgtDurCol);

% Duration scaling domain

durDom = opt.DurDomain;

if isempty(durDom) && hasDur
    dv = [edges.(srcDurCol); edges.(tgtDurCol)];
    dv = dv(isfinite(dv) & dv >= 0);
    if isempty(dv)
        durDom = [0 1];
    else
        durDom = [min(dv) max(dv)];
        if durDom(2) <= durDom(1)
            durDom = [durDom(1) durDom(1) + 1];
        end
    end
end

mapDurToRadius = @(d) local(d, durDom, opt.DurRadiusMin, opt.DurRadiusMax);

% --- Draw duration bubbles first (only for significant edges) ---

if hasDur && height(edges) > 0

    for r = 1: height(edges)

        i = edges.src(r);
        j = edges.tgt(r);

        if i == j || any([i j] < 1) || any([i j] > 7)
            continue;
        end

        p1 = [x(i) y(i)];
        p2 = [x(j) y(j)];
        v = p2 - p1;
        L = norm(v);

        if L <= eps
            continue;
        end

        u = v / L;

        dS = edges.(srcDurCol)(r);
        dT = edges.(tgtDurCol)(r);

        if ~isfinite(dS) || ~isfinite(dT)
            continue;
        end

        rS = mapDurToRadius(dS);
        rT = mapDurToRadius(dT);

        cS = p1 + u * (nodeR + rS);
        cT = p2 - u * (nodeR + rT);

        rectangle('Position', [cS(1) - rS, cS(2) - rS, 2 * rS, 2 * rS], 'Curvature', [1 1], 'FaceColor', opt.DurFaceColor, 'EdgeColor', opt.DurEdgeColor, 'LineWidth', opt.DurEdgeWidth);
        rectangle('Position', [cT(1) - rT, cT(2) - rT, 2 * rT, 2 * rT], 'Curvature', [1 1], 'FaceColor', opt.DurFaceColor, 'EdgeColor', opt.DurEdgeColor, 'LineWidth', opt.DurEdgeWidth);

    end
end

% --- Draw arrows (start/end constrained by duration bubbles when available) ---

for i = 1: 7
    for j = i + 1: 7

        has_ij = sig(i, j);
        has_ji = sig(j, i);

        if ~has_ij && ~has_ji
            continue;
        end

        p1 = [x(i) y(i)];
        p2 = [x(j) y(j)];

        if has_ij && ~ has_ji
            [s1, t1] = endpointsFromDur(edges, i, j, p1, p2, nodeR, hasDur, srcDurCol, tgtDurCol, mapDurToRadius);
            w1 = mapEffectToWidth(abs(log(irr(i, j))));
            col1 = iff(dirSign(i, j) <= 0, opt.BelowColor, opt.AboveColor);

            plot_filled_arrow([s1(1) t1(1)], [s1(2) t1(2)], w1, col1);

        end

        if has_ji && ~ has_ij
            [s2, t2] = endpointsFromDur(edges, j, i, p2, p1, nodeR, hasDur, srcDurCol, tgtDurCol, mapDurToRadius);
            w2 = mapEffectToWidth(abs(log(irr(j, i))));
            col2 = iff(dirSign(j, i) <= 0, opt.BelowColor, opt.AboveColor);

            plot_filled_arrow([s2(1) t2(1)], [s2(2) t2(2)], w2, col2);

        end

        if has_ij && has_ji

            [s1, t1] = endpointsFromDur(edges, i, j, p1, p2, nodeR, hasDur, srcDurCol, tgtDurCol, mapDurToRadius);
            [s2, t2] = endpointsFromDur(edges, j, i, p2, p1, nodeR, hasDur, srcDurCol, tgtDurCol, mapDurToRadius);

            w1 = mapEffectToWidth(abs(log(irr(i, j))));
            w2 = mapEffectToWidth(abs(log(irr(j, i))));

            col1 = iff(dirSign(i, j) <= 0, opt.BelowColor, opt.AboveColor);
            col2 = iff(dirSign(j, i) <= 0, opt.BelowColor, opt.AboveColor);

            drawBezierArrow_AbsGap(s1, t1, 0.8, col1, w1, 0, 0.2);
            drawBezierArrow_AbsGap(s2, t2, 0.8, col2, w2, 0, 0.2);

        end

    end
end

% Nodes on top

for k = 1: 7
    rectangle('Position', [x(k) - nodeR, y(k) - nodeR, 2 * nodeR, 2 * nodeR], 'Curvature', [1 1], 'EdgeColor', 'k', 'LineWidth', 2);
    text(x(k), y(k), labels{k}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 18);
end

hold off;

% ====================== helpers ======================

    function out = iff(cond, a, b)
        if cond
            out = a;
        else
            out = b;
        end
    end

    function r = local(d, dom, rmin, rmax)
        if ~isfinite(d)
            r = rmin;
            return;
        end

        if isempty(dom) || numel(dom) ~= 2 || dom(2) <= dom(1)
            r = 0.5 * (rmin + rmax);
            return;
        end

        t = (d - dom(1)) / (dom(2) - dom(1));
        t = min(1, max(0, t));
        r = rmin + t * (rmax - rmin);
    end

    function [pStart, pEnd] = endpointsFromDur(E, src, tgt, pSrc, pTgt, nodeRad, hasDurLocal, sCol, tCol, mapFn)

        v = pTgt - pSrc;
        L = norm(v);

        if L <= eps
            pStart = pSrc;
            pEnd = pTgt;
            return;
        end

        u = v / L;

        if ~hasDurLocal

            % Fallback: arrow touches node boundaries only

            pStart = pSrc + u * nodeRad;
            pEnd = pTgt - u * nodeRad;
            return;
        end

        idx = find(E.src == src & E.tgt == tgt, 1, 'first');

        if isempty(idx)
            pStart = pSrc + u * nodeRad;
            pEnd = pTgt - u * nodeRad;
            return;
        end

        dS = E.(sCol)(idx);
        dT = E.(tCol)(idx);

        if ~isfinite(dS) || ~isfinite(dT)
            pStart = pSrc + u * nodeRad;
            pEnd = pTgt - u * nodeRad;
            return;
        end

        rS = mapFn(dS);
        rT = mapFn(dT);

        cS = pSrc + u * (nodeRad + rS);
        cT = pTgt - u * (nodeRad + rT);

        pStart = cS + u * rS; % outer boundary of source duration bubble
        pEnd = cT - u * rT;% outer boundary of target duration bubble (facing source)
    end

    function drawBezierArrow_AbsGap(pStart, pEnd, curvFrac, col, lw, rnode, shaftGap)

        % Curved quadratic Bézier arrow from pStart -> pEnd.
        %
        % Behavior:
        % - Shaft starts at the source node boundary.
        % - Shaft ends an ABSOLUTE distance 'shaftGap' BEFORE the target node boundary.
        % - Arrowhead spans exactly that 'shaftGap' and its tip touches the target node.
        %
        % Note:
        % - shaftGap is in the same coordinate units as node positions.

        v = pEnd - pStart;
        L = norm(v);

        if L < eps
            return;
        end

        % Unit chord and normal

        u = v / L;
        n = [-u(2), u(1)];

        % Control point: midpoint + curvature bump (fraction of chord length)

        mid = (pStart + pEnd) / 2;
        ctrl = mid + curvFrac * (0.25 * L) * n;

        % Endpoint tangents for node trimming

        t0 = ctrl - pStart;
        t0 = t0 / max(eps, norm(t0));
        t1 = pEnd - ctrl;
        t1 = t1 / max(eps, norm(t1));

        % Trim endpoints to node boundaries

        a = pStart + t0 * (rnode * 1.05);
        b_tip = pEnd - t1 * (rnode * 1.05);

        % Sample the quadratic Bézier from a -> ctrl -> b_tip

        N = max(60, round(10 * L));
        t = linspace(0, 1, N).';
        A = (1 - t) .^ 2;
        B = 2 * (1 - t) .* t;
        C = t .^ 2;
        P = A .* a + B .* ctrl + C .* b_tip;
        dP = (-2 * (1 - t)) .* a + (2 - 4 * t) .* ctrl + (2 * t) .* b_tip;

        % Arc-length parameterization

        seg = sqrt(sum(diff(P, 1, 1) .^ 2, 2));
        s = [0; cumsum(seg)];
        totalLen = s(end);

        % Place shaft end at absolute distance 'shaftGap' before the tip

        gap = max(0, shaftGap);
        baseLen = max(0, totalLen - gap);

        % Interpolate base point (shaft end) at arc-length = baseLen

        idxBase = find(s >= baseLen, 1, 'first');

        if isempty(idxBase) || idxBase == 1
            basePt = P(1, :);
        else
            s0 = s(idxBase - 1);
            s1 = s(idxBase);
            lam = min(1, max(0, (baseLen - s0) / (s1 - s0 + eps)));
            basePt = (1 - lam) * P(idxBase - 1, :) + lam * P(idxBase, :);
        end

        % Draw shaft from start to basePt

        if baseLen > 0
            shaftX = [P(1: idxBase - 1, 1); basePt(1)];
            shaftY = [P(1: idxBase - 1, 2); basePt(2)];
            plot(shaftX, shaftY, 'Color', col, 'LineWidth', lw, 'Clipping', 'off');
        end

        % Arrowhead geometry: Make the arrowhead span the exact remaining gap to the tip

        tipTan = dP(end, :);

        if norm(tipTan) < eps
            tipTan = b_tip - basePt;
            if norm(tipTan) < eps
                tipTan = u;
            end
        end

        tipTan = tipTan / norm(tipTan);
        nvec = [-tipTan(2), tipTan(1)];

        % Use head length = shaftGap (clamped for reasonableness)

        head_len = min(max(gap, 0.02), 0.45);

        % Scale head width with both head_len and stroke width

        head_half_width = 0.5 * max(0.7 * head_len, 0.7 * lw * 0.02);

        % Place head base exactly at the shaft end (so it joins seamlessly)

        base_center = basePt;
        base_left = base_center + nvec * head_half_width * 1.4;
        base_right = base_center - nvec * head_half_width * 1.4;

        % Tip is at b_tip (on node boundary)

        patch([b_tip(1) base_left(1) base_right(1)], [b_tip(2) base_left(2) base_right(2)], col, 'EdgeColor', 'none', 'FaceAlpha', 1.0, 'Clipping', 'off');
    end


end