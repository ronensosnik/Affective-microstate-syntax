function edge_tbl = computeEdgeRatesTable(outputStats, windows, subj_id_map)
% computeEdgeRatesTable:  Tall table of per-edge rates per subject × condition × ERP.

rows = {};

grp_fields = fieldnames(outputStats);

cond_list = {'negative', 'neutral', 'positive'};

for g = 1: numel(grp_fields)
    grp_name = grp_fields{g};
    if any(strcmpi(grp_name, {'HC', 'Siblings'}))
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
                add_rows(rec, subjnum, group, mood, cond_name);
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
                add_rows(rec, subjnum, group, mood, cond_name);
            end
        end
    end
end

edge_tbl = vertcat(rows{ : });

    function add_rows(rec, subjnum, group, mood, cond_name)
        add_one(rec.MSClass, windows.N200, subjnum, group, mood, cond_name, 'N200');
        add_one(rec.MSClass, windows.P300, subjnum, group, mood, cond_name, 'P300');
        add_one(rec.MSClass, windows.LPP,  subjnum, group, mood, cond_name, 'LPP');
    end

    function add_one(MSClass, win, subjnum, group, mood, cond_name, ERPname)
        C = transitionCounts(MSClass, win);
        n_trials = size(MSClass, 2);
        wlen_s = (win(2) - win(1) + 1) ./ 1000;
        dur_sec = n_trials .* wlen_s;

        for i = 1: 7
            for j = 1: 7
                if i == j
                    continue
                end
                rate = C(i, j) ./ max(1e-9, dur_sec);
                T = table;
                T.subject = subjnum;
                T.group = string(group);
                T.mood = string(mood);
                T.condition = string(upperFirst(cond_name));
                T.ERP = string(ERPname);
                T.i = i;
                T.j = j;
                T.Edge = string([num2str(i) '->' num2str(j)]);
                T.edge_rate = rate;
                rows{end + 1, 1} = T;
            end
        end
    end

end
