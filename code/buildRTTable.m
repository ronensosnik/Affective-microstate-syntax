function rt_tbl = buildRTTable(France)

% buildRTTable:  Create a tidy table of RT summaries per subject × condition × ERP-stage placeholder
% Output columns:  subject, group, mood, condition, median_rt_ms, log_median_rt, n_used, age, gender, MADRS, YMRS, GAF

subj_vec = [];
group_vec = strings(0, 1);
mood_vec = strings(0, 1);
cond_vec = strings(0, 1);
median_rt_ms_vec = [];
log_median_rt_vec = [];
n_used_vec = [];
age_vec = [];
gender_vec = strings(0, 1);
MADRS_vec = [];
YMRS_vec = [];
GAF_vec = [];

process_group_mood(France, 'BP_I', 'Depressed');
process_group_mood(France, 'BP_I', 'Euthymic');
process_group_mood(France, 'BP_II', 'Depressed');
process_group_mood(France, 'BP_II', 'Euthymic');
process_group_nomood(France, 'HC');
process_group_nomood(France, 'Siblings');

rt_tbl = table(subj_vec, group_vec, mood_vec, cond_vec, median_rt_ms_vec, log_median_rt_vec, n_used_vec, age_vec, gender_vec, MADRS_vec, YMRS_vec, GAF_vec, 'VariableNames', {'subject', 'group', 'mood', 'condition', 'median_rt_ms', 'log_median_rt', 'n_used', 'age', 'gender', 'MADRS', 'YMRS', 'GAF'});

    function process_group_mood(F, group_name, mood_name)

        % process_group_mood:  Iterate subjects for groups with mood subfolders

        if ~isfield(F, group_name)
            return
        end

        if ~isfield(F.(group_name), mood_name)
            return
        end

        cells = F.(group_name).(mood_name);

        for k = 1: numel(cells)

            s = cells{k};

            if isempty(s)
                continue
            end

            add_subject_rows(s, group_name, mood_name);
        end
    end

    function process_group_nomood(F, group_name)

        % process_group_nomood:  Iterate subjects for groups without mood subfolders

        if ~isfield(F, group_name)
            return
        end

        cells = F.(group_name);

        for k = 1: numel(cells)
            s = cells{k};

            if isempty(s)
                continue
            end

            add_subject_rows(s, group_name, 'NA');
        end
    end

    function add_subject_rows(s, group_name, mood_name)

        % add_subject_rows:  Append one row per condition with RT summary, ignoring NaNs

        add_one_condition(s, group_name, mood_name, 'Neutral', s.Neutral);
        add_one_condition(s, group_name, mood_name, 'Negative', s.Negative);
        add_one_condition(s, group_name, mood_name, 'Positive', s.Positive);
    end

    function add_one_condition(s, group_name, mood_name, cond_name, cond_struct)

        % add_one_condition:  Compute median RT and add to vectors

        rt = cond_struct.Correct_resp_RT;
        rt = rt(isfinite(rt));

        if isempty(rt)
            return
        end

        subj_vec(end + 1, 1) = s.number;
        group_vec(end + 1, 1) = string(group_name);
        mood_vec(end + 1, 1) = string(mood_name);
        cond_vec(end + 1, 1) = string(cond_name);
        med_rt = median(rt);
        median_rt_ms_vec(end + 1, 1) = med_rt;
        log_median_rt_vec(end + 1, 1) = log(med_rt);
        n_used_vec(end + 1, 1) = numel(rt);
        age_vec(end + 1, 1) = s.age;
        gender_vec(end + 1, 1) = string(s.gender);

        if isfield(s, 'MADRS')
            MADRS_vec(end + 1, 1) = s.MADRS;
        else
            MADRS_vec(end + 1, 1) = NaN;
        end

        if isfield(s, 'YMRS')
            YMRS_vec(end + 1, 1) = s.YMRS;
        else
            YMRS_vec(end + 1, 1) = NaN;
        end

        if isfield(s, 'GAF')
            GAF_vec(end + 1, 1) = s.GAF;
        else
            GAF_vec(end + 1, 1) = NaN;
        end

    end
end