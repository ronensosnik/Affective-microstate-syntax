function drawERPGridFigures(TransStats, erp, grp, Conds, varargin)

p = inputParser;

% ====================================================================%

% Figure name

addParameter(p, 'FigureName', sprintf('Significant transitions — %s | %s', erp, grp));

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

nCols = numel(Conds);
figure('Name', opt.FigureName, 'Color', 'w', 'Units', 'normalized', 'Position', [0.0500 0.6022 0.2750 0.1489]);
t = tiledlayout(1, nCols, 'Padding', 'none', 'TileSpacing', 'none');

for c = 1: numel(Conds)

    condName = Conds{c};
    nexttile(t, c);

    if isfield(TransStats.(erp).(grp), condName)
        E = TransStats.(erp).(grp).(condName).edges;
    else
        E = table([], [], [], [], [], [], strings(0, 1), [], [], 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'ObsOverExp', 'Direction', 'source_dur', 'target_dur'});
    end

    plotCircularMicrostateTransitions(E, 'NodeLabels', opt.NodeLabels, 'NodeRadius', opt.NodeRadius, 'GreyMesh', 'true', 'GreyLineWidth', opt.GreyLineWidth, 'GreyColor', opt.GreyColor, 'WidthDomain', opt.WidthDomain, 'MinArrowWidth', opt.MinArrowWidth, 'MaxArrowWidth', opt.MaxArrowWidth, 'AboveColor', opt.AboveColor, 'BelowColor', opt.BelowColor, 'DurDomain', opt.DurDomain, 'DurRadiusMin', opt.DurRadiusMin, 'DurRadiusMax', opt.DurRadiusMax, 'DurFaceColor', opt.DurFaceColor, 'DurEdgeColor', opt.DurEdgeColor, 'DurEdgeWidth', opt.DurEdgeWidth);
end

end