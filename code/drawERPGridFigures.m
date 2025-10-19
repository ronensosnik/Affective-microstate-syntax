function drawERPGridFigures(TransStats, erp, grp, Conds, MicrostateLabels, varargin) % One figure per ERP×Group with 3 large circular plots (one per condition). 

p = inputParser; 

 addParameter(p, 'FigureName', sprintf('Significant transitions — %s | %s', erp, grp)); 
 addParameter(p, 'NodeRadius', []); 
 addParameter(p, 'ArrowOffset', []); 
 addParameter(p, 'WidthDomain', []); 
 addParameter(p, 'CurveMagnitude', []); 
 addParameter(p, 'MinWidth', []); 
 addParameter(p, 'MaxWidth', []); 
 addParameter(p, 'GreyLineWidth', []); 
 addParameter(p, 'GreyColor', []);
 addParameter(p, 'MicrostateLabels', []);
 addParameter(p, 'ShaftGap', []);
 addParameter(p, 'AboveColor', []);
 addParameter(p, 'BelowColor', []);

 addParameter(p, 'GreyMesh', true); 

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
         E = table([], [], [], [], [], [], strings(0, 1), 'VariableNames', {'src', 'tgt', 'IRR', 'IRR_CIlo', 'IRR_CIhi', 'ObsOverExp', 'Direction'});
     end 
     
     plotCircularMicrostateTransitions(E, 'Labels', MicrostateLabels, 'NodeRadius', opt.NodeRadius, 'ArrowOffset', opt.ArrowOffset, 'GreyMesh', opt.GreyMesh, 'GreyLineWidth', opt.GreyLineWidth, 'GreyColor', opt.GreyColor, 'WidthDomain', opt.WidthDomain, 'MinWidth', opt.MinWidth, 'MaxWidth', opt.MaxWidth, 'CurveMagnitude', opt.CurveMagnitude, 'ShaftGap', opt.ShaftGap, 'AboveColor', opt.AboveColor, 'BelowColor', opt.BelowColor); 
 end
end