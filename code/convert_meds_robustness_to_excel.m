function [xlsxPath, sheetsWritten] = convert_meds_robustness_to_excel(rootDir, outFile)
% convert_meds_robustness_to_excel
% Auto-discovers the six "*_Baseline_vs_meds[ _poly].mat" files, loads their
% .detail and .summary tables, and writes them to "Metrics_Meds_Robustness.xlsx".
%
% This reproduces the workbook format:
%   <Measure>_<classes|poly>_detail
%   <Measure>_<classes|poly>_summary
% for Measure in {"Coverage","Duration","Occurrence"}.
%
% USAGE
%   convert_meds_robustness_to_excel;                 % search from pwd
%   convert_meds_robustness_to_excel('C:\path\to\dir'); % search from given dir
%   convert_meds_robustness_to_excel('C:\dir','C:\out\Metrics_Meds_Robustness.xlsx');
%
% OUTPUTS
%   xlsxPath      -> full path to the written workbook
%   sheetsWritten -> struct with sheet names and sizes
%
% NOTES
% - The .mat files are expected to contain a single variable whose value is a
%   struct with fields 'detail' and 'summary' (tables). This is exactly what
%   summarize_metric_meds_robustness saved.
% - If your MATLAB is pre-R2016b (no "**" in dir), the function falls back to
%   genpath-based recursion.

if nargin < 1 || isempty(rootDir), rootDir = pwd; end
if nargin < 2 || isempty(outFile), outFile = fullfile(rootDir,'Metrics_Meds_Robustness.xlsx'); end
xlsxPath = outFile;

measures = {'Coverage','Duration','Occurrence'};
specs    = {'classes','poly'};  % 'classes' => file without "_poly", 'poly' => file with "_poly"

% Make sure output dir exists
outDir = fileparts(outFile);
if ~isempty(outDir) && ~exist(outDir,'dir')
    mkdir(outDir);
end

% If the workbook already exists, start fresh so old sheets don't linger.
if exist(outFile,'file')
    try
        delete(outFile);
    catch ME
        warning('Could not delete existing %s: %s', outFile, ME.message);
    end
end

sheetsWritten = struct;
fprintf('[convert_meds_robustness_to_excel] Writing: %s\n', outFile);

for m = 1:numel(measures)
    measure = measures{m};

    for s = 1:numel(specs)
        spec = specs{s};
        % Build expected filename
        if strcmpi(spec,'classes')
            baseName = sprintf('%s_Baseline_vs_meds.mat', measure);
        else
            baseName = sprintf('%s_Baseline_vs_meds_poly.mat', measure);
        end

        % Find the file (recursive)
        matPath = findFirstFile(rootDir, baseName);
        if isempty(matPath)
            error('Required file not found: %s (searched under "%s")', baseName, rootDir);
        end

        % Load the single variable inside
        S = load(matPath);
        varNames = fieldnames(S);
        if isempty(varNames)
            error('No variables found in %s', matPath);
        end
        data = S.(varNames{1});

        % Expect fields 'detail' and 'summary'
        if ~isfield(data,'detail') || ~isfield(data,'summary')
            error('File %s does not contain fields "detail" and "summary".', matPath);
        end

        detailTbl  = ensureTable(data.detail);
        summaryTbl = ensureTable(data.summary);

        % Compose sheet names to match prior export
        shDetail  = sprintf('%s_%s_detail',  measure, spec);
        shSummary = sprintf('%s_%s_summary', measure, spec);

        % Write the two sheets
        writetable(detailTbl,  outFile, 'Sheet', shDetail);
        writetable(summaryTbl, outFile, 'Sheet', shSummary);

        % Track sizes for console/reporting
        sheetsWritten.(shDetail)  = [height(detailTbl)  width(detailTbl)];
        sheetsWritten.(shSummary) = [height(summaryTbl) width(summaryTbl)];

        fprintf('  - %s: wrote %-28s (%d x %d)\n', measure, shDetail,  height(detailTbl),  width(detailTbl));
        fprintf('  - %s: wrote %-28s (%d x %d)\n', measure, shSummary, height(summaryTbl), width(summaryTbl));
    end
end

% Add a tiny README sheet at the end (optional)
try
    readme = {
        'GeneratedBy', 'convert_meds_robustness_to_excel.m'
        'GeneratedOn', datestr(now, 'yyyy-mm-dd HH:MM:SS')
        'SourceDir',   char(rootDir)
        'Notes',       'Each measure has "classes" (per-class) and "poly" (poly/pooled) specs, each with detail and summary sheets.'
        };
    writecell(readme, outFile, 'Sheet', 'README', 'Range', 'A1');
catch ME
    warning('README sheet not written: %s', ME.message);
end

fprintf('[convert_meds_robustness_to_excel] Done.\n');

end % function main


% ---------- helpers ----------

function T = ensureTable(x)
% Convert dataset/struct/whatever into a table if needed.
    if istable(x)
        T = x;
        return;
    end

    % Support older Statistics Toolbox "dataset"
    if isa(x,'dataset')
        try
            T = dataset2table(x);
            return;
        catch
            % fall through
        end
    end

    % Struct array -> table
    if isstruct(x)
        try
            T = struct2table(x);
            return;
        catch
            % If it's a scalar struct of tables, concatenate fields horizontally
            f = fieldnames(x);
            if isscalar(x) && all(cellfun(@(fn) istable(x.(fn)), f))
                % Concatenate side-by-side with <field>_ prefixes if needed
                T = x.(f{1});
                for i = 2:numel(f)
                    Ti = x.(f{i});
                    % Make unique names with prefix
                    Ti.Properties.VariableNames = matlab.lang.makeUniqueStrings( ...
                        strcat(f{i}, "_", string(Ti.Properties.VariableNames)), ...
                        [T.Properties.VariableNames]);
                    T = [T Ti]; %#ok<AGROW>
                end
                return;
            end
        end
    end

    % Final fallback: make a single-cell table
    warning('Object of type %s was not a table; writing as a single-cell table.', class(x));
    T = table({x}, 'VariableNames', {'Value'});
end

function fpath = findFirstFile(rootDir, fileName)
% Recursive search for fileName under rootDir (first match).
    fpath = '';
    % Try fast modern recursive dir (R2016b+)
    try
        D = dir(fullfile(rootDir, '**', fileName));
    catch
        % Fallback for older MATLAB: manual recursion via genpath
        D = [];
        p = genpath(rootDir);
        parts = strsplit(p, pathsep);
        for i = 1:numel(parts)
            if isempty(parts{i}), continue; end
            Di = dir(fullfile(parts{i}, fileName));
            if ~isempty(Di)
                D = Di(1); %#ok<AGROW>
                break;
            end
        end
        if ~isempty(D)
            fpath = fullfile(D.folder, D.name);
            return;
        end
    end

    if ~isempty(D)
        fpath = fullfile(D(1).folder, D(1).name);
    end
end
