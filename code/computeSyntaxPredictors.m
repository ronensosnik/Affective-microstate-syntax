function syntax_tbl = computeSyntaxPredictors(outputStats, anchors, windows, subj_id_map)

% computeSyntaxPredictors:  Per subject × condition × ERP motif predictors.
%
% Returns columns:
% subject, group, mood, condition, ERP, anchor1_rate, anchor2_rate, backbone_rate,
% anchor1_delta, anchor2_delta, backbone_delta, n_trials, window_sec

rows = {};

grp_fields = fieldnames(outputStats);

for g = 1: numel(grp_fields)

    grp_name = grp_fields{g};
    cond_list = {'negative', 'neutral', 'positive'};

    if strcmpi(grp_name, 'HC') || strcmpi(grp_name, 'Siblings')

        group = grp_name;
        mood = 'NA';
        ids_vec = subj_id_map.(grp_name);

        for c = 1: numel(cond_list)
            cond_name = cond_list{c};
            subj_arr = outputStats.(grp_name).(cond_name);

            for s = 1: numel(subj_arr)
                rec = subj_arr(s);
                subjnum = NaN;

                if s <= numel(ids_vec)
                    subjnum = ids_vec(s);
                end

                add_rows(rec, subjnum, group, mood, cond_name, anchors, windows);
            end
        end
    else

        if startsWith(grp_name, 'BP_I_')
            group = 'BP_I';
            mood = erase(grp_name, 'BP_I_');

        elseif startsWith(grp_name, 'BP_II_')
            group = 'BP_II';
            mood = erase(grp_name, 'BP_II_');

        else
            continue
        end

        ids_vec = subj_id_map.(grp_name);

        for c = 1: numel(cond_list)

            cond_name = cond_list{c};
            subj_arr = outputStats.(grp_name).(cond_name);

            for s = 1: numel(subj_arr)

                rec = subj_arr(s);
                subjnum = NaN;

                if s <= numel(ids_vec)
                    subjnum = ids_vec(s);
                end

                add_rows(rec, subjnum, group, mood, cond_name, anchors, windows);
            end
        end
    end
end

syntax_tbl = vertcat(rows{ : });

    function add_rows(rec, subjnum, group, mood, cond_name, anchors, windows)

        rowN = compute_row(rec.MSClass, windows.N200, anchors.N200, [], subjnum, group, mood, cond_name, 'N200');
        rows{end + 1, 1} = rowN;

        rowP = compute_row(rec.MSClass, windows.P300, anchors.P300, [], subjnum, group, mood, cond_name, 'P300');
        rows{end + 1, 1} = rowP;

        rowL = compute_row(rec.MSClass, windows.LPP, [], anchors.LPP, subjnum, group, mood, cond_name, 'LPP');
        rows{end + 1, 1} = rowL;
    end

    function T = compute_row(MSClass, win, anchor_pair, backbone_pair, subjnum, group, mood, cond_name, ERPname)

        C = transitionCounts(MSClass, win);
        occ = occupancyCounts(MSClass, win);

        Ttot = sum(C(:));
        p = occ ./ max(1, sum(occ));

        n_trials = size(MSClass, 2);
        wlen_s = (win(2) - win(1) + 1) ./ 1000;
        dur_sec = n_trials .* wlen_s;

        anchor1_rate = NaN;
        anchor2_rate = NaN;
        backbone_rate = NaN;

        anchor1_delta = NaN;
        anchor2_delta = NaN;
        backbone_delta = NaN;

        if ~isempty(anchor_pair)
            a1 = anchor_pair(1);
            a2 = anchor_pair(2);

            sumObs1 = sum(C(a1, :)) - C(a1, a1);
            sumObs2 = sum(C(a2, :)) - C(a2, a2);

            anchor1_rate = sumObs1 ./ max(1e-9, dur_sec);
            anchor2_rate = sumObs2 ./ max(1e-9, dur_sec);

            sumExp1 = Ttot .* p(a1) .* (1 - p(a1));
            sumExp2 = Ttot .* p(a2) .* (1 - p(a2));

            if sumExp1 > 0
                anchor1_delta = (sumObs1 - sumExp1) ./ sumExp1;
            else
                anchor1_delta = NaN;
            end

            if sumExp2 > 0
                anchor2_delta = (sumObs2 - sumExp2) ./ sumExp2;
            else
                anchor2_delta = NaN;
            end
        end

        if ~isempty(backbone_pair)
            i = backbone_pair(1);
            j = backbone_pair(2);

            obs_ij = C(i, j);
            obs_ji = C(j, i);

            backbone_rate = ((obs_ij ./ max(1e-9, dur_sec)) + (obs_ji ./ max(1e-9, dur_sec))) ./ 2;

            sumExp = 2 .* Ttot .* p(i) .* p(j);

            if sumExp > 0
                backbone_delta = ((obs_ij + obs_ji) - sumExp) ./ sumExp;
            else
                backbone_delta = NaN;
            end
        end

        T = table;
        T.subject = subjnum;
        T.group = string(group);
        T.mood = string(mood);
        T.condition = string(upperFirst(cond_name));
        T.ERP = string(ERPname);
        T.anchor1_rate = anchor1_rate;
        T.anchor2_rate = anchor2_rate;
        T.backbone_rate = backbone_rate;
        T.anchor1_delta = anchor1_delta;
        T.anchor2_delta = anchor2_delta;
        T.backbone_delta = backbone_delta;
        T.n_trials = n_trials;
        T.window_sec = dur_sec;
    end

end