function drawERPGridFigures(TransStats, erp, grp, Conds, MicrostateLabels, varargin) % One figure per ERP×Group with 3 large circular plots (one per condition). 

%% drawERPGridFigures.m
% Pipeline to back-fit 7-class Custo2017 templates, segment N200/P300/LPP windows,
% and compute duration, coverage, occurrence, and directed transition counts.
% 
% Inputs (via config/config.m):
%   data_deriv : folder with preprocessed/epoched .set files organized by Group/Subject/Condition
%   results_dir: root output
% Dependencies:
%   EEGLAB, MicrostateLab (Custo2017), Statistics & ML Toolbox.
% Outputs:
%   results/stats/outputStats.mat      % per-subject microstate sequences & GFP
%   results/stats/OUT_*.mat            % OUT_Coverage / OUT_Occurrence / OUT_Duration
%   results/figures/*.png              % metrics and transition plots
% Notes:
%   - Uses continuous back-fitting (not GFP-peaks only); smoothing (lambda, b) adjustable.
%   - Windows: N200 180–300 ms; P300 300–500 ms; LPP 500–1000 ms.
%   - Replace any absolute paths with variables from config.m.

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
