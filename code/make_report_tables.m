function REP = make_report_tables(OUT_measure)
% Build publication-ready tables from OUT = lme_posthoc_21FDR_transformsBT(...)

ERPs   = OUT_measure.info.ERPs;
Micros = OUT_measure.info.Microstates;
measure= OUT_measure.info.measure;

Omni = table();   % omnibus (interaction, group, condition) + FDR + model used
AgeT = table();   % Age effect (F-test on chosen model) + df + p
EMMs = table();   % back-transformed estimated marginal means (+95% CI, units, Ns)
Pairs= table();   % pairwise comparisons with back-transformed effect sizes + FDR

for e = 1:numel(ERPs)
  for m = 1:numel(Micros)
    R = OUT_measure.(measure){e,m};
    if isempty(R), continue; end

    % ---------- Omnibus row ----------
    omniRow = table( ...
      string(measure), string(ERPs{e}), string(Micros{m}), string(R.model_used), ...
      R.tests.interaction.p, R.tests.interaction.pFDR, ...
      R.tests.group.p,       R.tests.group.pFDR, ...
      R.tests.condition.p,   R.tests.condition.pFDR, ...
      'VariableNames', {'Measure','ERP','Microstate','ModelUsed', ...
                        'p_Interaction','pFDR_Interaction', ...
                        'p_Group','pFDR_Group', 'p_Condition','pFDR_Condition'});
    Omni = [Omni; omniRow];

    % ---------- Age effect on chosen model ----------
    % pick the model used for inference
    switch R.model_used
      case 'full_with_interaction'; M = R.models.full_with_interaction;
      otherwise;                    M = R.models.no_interaction_REML;
    end
    A = anova(M,'DFMethod','Satterthwaite');
    if any(strcmp(A.Term,'Age_c'))
      rowA = A(strcmp(A.Term,'Age_c'),:);
      ageRow = table(string(measure), string(ERPs{e}), string(Micros{m}), ...
        rowA.DF1, rowA.DF2, rowA.FStat, rowA.pValue, ...
        'VariableNames',{'Measure','ERP','Microstate','DF1','DF2','F_Age','p_Age'});
      AgeT = [AgeT; ageRow];
    else
      ageRow = table(string(measure), string(ERPs{e}), string(Micros{m}), ...
        NaN,NaN,NaN,NaN, 'VariableNames',{'Measure','ERP','Microstate','DF1','DF2','F_Age','p_Age'});
      AgeT = [AgeT; ageRow];
    end

    % ---------- EMMs (back-transformed) ----------
    % Pull whichever EMMs exist in R (depends on significance branch)
    emmTables = {};
    if isfield(R,'emm')
      fns = fieldnames(R.emm);
      for k=1:numel(fns)
        T = R.emm.(fns{k});
        if ~isempty(T)
          T2 = T;  % add labels
          T2.Measure   = repmat(string(measure),height(T2),1);
          T2.ERP       = repmat(string(ERPs{e}),height(T2),1);
          T2.Microstate= repmat(string(Micros{m}),height(T2),1);
          T2.Source    = repmat(string(fns{k}),height(T2),1); % e.g., Group_by_Condition
          % Add Ns per cell (unique SubjectID per level/By-level)
          if ismember('By', lower(T.Properties.VariableNames))
            % not needed; margins returns 'By' as separate columns
          end
          % compute Ns from R.dataTable:
          if strcmpi(fns{k},'Group_by_Condition')
            G = categories(R.dataTable.Group); C = categories(R.dataTable.Condition);
            Ncol = nan(height(T2),1);
            for r=1:height(T2)
              g = T2.Group(r); c = T2.Condition(r);
              mask = (R.dataTable.Group==g) & (R.dataTable.Condition==c);
              Ncol(r) = numel(unique(R.dataTable.SubjectID(mask)));
            end
            T2.N_subj = Ncol;
          elseif strcmpi(fns{k},'Condition_by_Group')
            G = categories(R.dataTable.Group); C = categories(R.dataTable.Condition);
            Ncol = nan(height(T2),1);
            for r=1:height(T2)
              g = T2.Group(r); c = T2.Condition(r);
              mask = (R.dataTable.Group==g) & (R.dataTable.Condition==c);
              Ncol(r) = numel(unique(R.dataTable.SubjectID(mask)));
            end
            T2.N_subj = Ncol;
          elseif strcmpi(fns{k},'Group')
            G = categories(R.dataTable.Group);
            Ncol = nan(height(T2),1);
            for r=1:height(T2)
              g = T2.Group(r);
              mask = (R.dataTable.Group==g);
              Ncol(r) = numel(unique(R.dataTable.SubjectID(mask)));
            end
            T2.N_subj = Ncol;
          elseif strcmpi(fns{k},'Condition')
            C = categories(R.dataTable.Condition);
            Ncol = nan(height(T2),1);
            for r=1:height(T2)
              c = T2.Condition(r);
              mask = (R.dataTable.Condition==c);
              Ncol(r) = numel(unique(R.dataTable.SubjectID(mask)));
            end
            T2.N_subj = Ncol;
          end
          % keep only readable columns
          keep = intersect({'Measure','ERP','Microstate','Source', ...
                            'Group','Condition', ...
                            'EMM_BT','EMM_BT_Low','EMM_BT_High','Units','N_subj'}, ...
                            T2.Properties.VariableNames,'stable');
          emmTables{end+1} = T2(:,keep);
        end
      end
    end
    if ~isempty(emmTables), EMMs = [EMMs; vertcat(emmTables{:})]; end

    % ---------- Pairwise (back-transformed + FDR) ----------
    if isfield(R,'pairwise')
      fns = fieldnames(R.pairwise);
      for k=1:numel(fns)
        P = R.pairwise.(fns{k});
        if ~isempty(P)
          P2 = P;
          P2.Measure    = repmat(string(measure),height(P2),1);
          P2.ERP        = repmat(string(ERPs{e}),height(P2),1);
          P2.Microstate = repmat(string(Micros{m}),height(P2),1);
          P2.Source     = repmat(string(fns{k}),height(P2),1); % e.g., Group_by_Condition, Group, Condition

          % Normalize pair labels (MATLAB names depend on factor)
          vnames = P2.Properties.VariableNames;
          g1 = vnames(startsWith(vnames,'Group_'));  % e.g., {'Group_1','Group_2'}
          c1 = vnames(startsWith(vnames,'Condition_'));
          if numel(g1)==2, P2.Level1 = P2.(g1{1}); P2.Level2 = P2.(g1{2}); end
          if numel(c1)==2, P2.Level1 = P2.(c1{1}); P2.Level2 = P2.(c1{2}); end

          % Keep readable columns depending on measure
          baseCols = {'Measure','ERP','Microstate','Source','Level1','Level2'};
          if strcmpi(measure,'Duration')
              effCols = {'Delta','Delta_Low','Delta_High'};
              effName = {'Effect','Effect_low','Effect_high'};
              units   = "ms";
          elseif strcmpi(measure,'Coverage')
              effCols = {'OR','OR_Low','OR_High'};
              effName = {'OR','OR_low','OR_high'};
              units   = "odds ratio";
          else
              effCols = {'RR','RR_Low','RR_High'};
              effName = {'RR','RR_low','RR_high'};
              units   = "rate ratio";
          end
          hasEff = all(ismember(effCols, P2.Properties.VariableNames));
          keep = [baseCols, {'Estimate','SE','DF','pValue','pFDR'}];
          if hasEff, keep = [keep, effCols]; end
          P2 = P2(:, intersect(keep, P2.Properties.VariableNames,'stable'));
          if hasEff
              P2.Properties.VariableNames(end-2:end) = effName;
              if ~ismember('Units', P2.Properties.VariableNames)
                  P2.Units = repmat(units,height(P2),1);
              end
          end
          Pairs = [Pairs; P2];
        end
      end
    end

  end
end

% Nice ordering
if ~isempty(Omni)
  Omni = movevars(Omni, {'Measure','ERP','Microstate','ModelUsed'}, 'Before',1);
end
if ~isempty(AgeT)
  AgeT = movevars(AgeT, {'Measure','ERP','Microstate'}, 'Before',1);
end
if ~isempty(EMMs)
  EMMs = movevars(EMMs, {'Measure','ERP','Microstate','Source'}, 'Before',1);
end
if ~isempty(Pairs)
  Pairs = movevars(Pairs, {'Measure','ERP','Microstate','Source','Level1','Level2'}, 'Before',1);
end

REP = struct('Omnibus',Omni,'Age',AgeT,'EMMs',EMMs,'Pairs',Pairs);
end
