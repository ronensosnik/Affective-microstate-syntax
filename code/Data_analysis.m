close all; clear all; clc;

check_Dataset = 0;

Compute_performance_indices_each_subject = 0;
Generate_performance_table_all_subjects = 1;
Analyze_data_psychophysics = 0;

Group_paper = {'BP_I_Depressed', 'BP_I_Euthymic', 'BP_II_Depressed', 'BP_II_Euthymic', 'HC', 'Siblings'};
Cond = {'Neutral', 'Negative', 'Positive'};

France_subj_BP_I_Depressed = [13 14 26 40 45 49 54 59 64];
France_subj_BP_I_Euthymic = [16 19 28 32 58 60 68 70 82];
France_subj_BP_II_Depressed = [24 29 35 36 37 50 53 55 72 77 83 84];
France_subj_BP_II_Euthymic = [7 10 15 33 48 52 56 75 76 80 81 88];
France_subj_HC = [20 23 27 30 39 41 42 47 51 62 71 73];
France_subj_Siblings = [6 9 17 22 38 65 66 67 69 79 86 87];

if Compute_performance_indices_each_subject

    %% Read .set data

    % Here, we generate a databasec contatining the response time of all groups in all conditions. We do not yet take off trials in which
    % the RT was too early (RT < 150ms) or that had no response in order to have the full dataset for later generating a table of excluded trials

    for group = 1: length(Group_paper)
        try

            cd (['D:\Papers\2025\In preparation\XXX (task microstates)\Data\' Group_paper{group}]);

            a = dir;

            for n = 1: length(a)

                nn_Neutral = 0; nn_Negative = 0; nn_Positive = 0;

                if a(n).isdir && ~strcmp(a(n).name, '.') && ~strcmp(a(n).name, '..')

                    try

                        cd (['D:\Papers\2025\In preparation\XXX (task microstates)\Data\' Group_paper{group} '\' a(n).name]);
                        dd = dir('*CRF.mat');
                        clear Subject; load (dd.name);

                        cd (['D:\Papers\2025\In preparation\XXX (task microstates)\Data\' Group_paper{group} '\' a(n).name]);
                        inputname = [a(n).name '_task_processed.set']; inputpath = pwd;

                        clear EEG command; [EEG, command] = pop_loadset(inputname, inputpath, []);

                        Subject.Neutral.Correct_RT = []; Subject.Neutral.Incorrect_RT = []; Subject.Neutral.Correct_inCorrect_index = []; Subject.Neutral.NoResponse_index = [];
                        Subject.Negative.Correct_RT = []; Subject.Negative.Incorrect_RT = []; Subject.Negative.Correct_inCorrect_index = []; Subject.Negative.NoResponse_index = [];
                        Subject.Positive.Correct_RT = []; Subject.Positive.Incorrect_RT = []; Subject.Positive.Correct_inCorrect_index = []; Subject.Positive.NoResponse_index = [];

                        for trial = 1: size(EEG.event, 2) - 1

                            if strcmp(EEG.event(trial).type,'100') || strcmp(EEG.event(trial).type,'101') || strcmp(EEG.event(trial).type, 'S100') || strcmp(EEG.event(trial).type, 'S101')

                                nn_Neutral = nn_Neutral + 1;

                                if (strcmp(EEG.event(trial + 1).type, '250') || strcmp(EEG.event(trial + 1).type, 'S250'))
                                    if round(EEG.event(trial + 1).latency - EEG.event(trial).latency) >= 150 && round(EEG.event(trial + 1).latency - EEG.event(trial).latency) <= 1995

                                        % There was a response, it was in the allowed time window, it was correct

                                        Subject.Neutral.TooEarly(nn_Neutral) = 0/0;

                                        Subject.Neutral.Correct_inCorrect_index(nn_Neutral) = 1;
                                        Subject.Neutral.NoResponse_index(nn_Neutral) = 0/0;

                                        Subject.Neutral.Correct_RT(nn_Neutral) = round(EEG.event(trial + 1).latency - EEG.event(trial).latency);
                                        Subject.Neutral.Incorrect_RT(nn_Neutral) = 0/0;

                                    elseif round(EEG.event(trial + 1).latency - EEG.event(trial).latency) < 150

                                        % There was a response, it was too early, it was correct

                                        Subject.Neutral.TooEarly(nn_Neutral) = 1;

                                        Subject.Neutral.Correct_inCorrect_index(nn_Neutral) = 0/0;
                                        Subject.Neutral.NoResponse_index(nn_Neutral) = 0/0;

                                        Subject.Neutral.Correct_RT(nn_Neutral) = 0/0;
                                        Subject.Neutral.Incorrect_RT(nn_Neutral) = 0/0;

                                    elseif (round(EEG.event(trial + 1).latency - EEG.event(trial).latency) > 1995)

                                        % There was no response

                                        Subject.Neutral.TooEarly(nn_Neutral) = 0/0;

                                        Subject.Neutral.Correct_inCorrect_index(nn_Neutral) = 0/0;
                                        Subject.Neutral.NoResponse_index(nn_Neutral) = 1;

                                        Subject.Neutral.Correct_RT(nn_Neutral) = 0/0;
                                        Subject.Neutral.Incorrect_RT(nn_Neutral) = 0/0;

                                    end

                                elseif strcmp(EEG.event(trial + 1).type, '251') || strcmp(EEG.event(trial + 1).type, 'S251')

                                    if round(EEG.event(trial + 1).latency - EEG.event(trial).latency) >= 150 && round(EEG.event(trial + 1).latency - EEG.event(trial).latency) <= 1995

                                        % There was a response, it was in the allowed time window, it was incorrect

                                        Subject.Neutral.TooEarly(nn_Neutral) = 0/0;

                                        Subject.Neutral.Correct_inCorrect_index(nn_Neutral) = 0;
                                        Subject.Neutral.NoResponse_index(nn_Neutral) = 0/0;

                                        Subject.Neutral.Correct_RT(nn_Neutral) = round(EEG.event(trial + 1).latency - EEG.event(trial).latency);
                                        Subject.Neutral.Incorrect_RT(nn_Neutral) = 0/0;

                                    elseif round(EEG.event(trial + 1).latency - EEG.event(trial).latency) < 150

                                        % There was a response, it was too early, it was incorrect

                                        Subject.Neutral.TooEarly(nn_Neutral) = 1;

                                        Subject.Neutral.Correct_inCorrect_index(nn_Neutral) = 0/0;
                                        Subject.Neutral.NoResponse_index(nn_Neutral) = 0/0;

                                        Subject.Neutral.Correct_RT(nn_Neutral) = 0/0;
                                        Subject.Neutral.Incorrect_RT(nn_Neutral) = 0/0;

                                    elseif (round(EEG.event(trial + 1).latency - EEG.event(trial).latency) > 1995)

                                        % There was no response

                                        Subject.Neutral.TooEarly(nn_Neutral) = 0/0;

                                        Subject.Neutral.Correct_inCorrect_index(nn_Neutral) = 0/0;
                                        Subject.Neutral.NoResponse_index(nn_Neutral) = 1;

                                        Subject.Neutral.Correct_RT(nn_Neutral) = 0/0;
                                        Subject.Neutral.Incorrect_RT(nn_Neutral) = 0/0;

                                    end
                                end

                            elseif strcmp(EEG.event(trial).type, '150') || strcmp(EEG.event(trial).type, '151') || strcmp(EEG.event(trial).type, 'S150') || strcmp(EEG.event(trial).type, 'S151')

                                nn_Negative = nn_Negative + 1;

                                if (strcmp(EEG.event(trial + 1).type, '250') || strcmp(EEG.event(trial + 1).type, 'S250'))
                                    if round(EEG.event(trial + 1).latency - EEG.event(trial).latency) >= 150 && round(EEG.event(trial + 1).latency - EEG.event(trial).latency) <= 1995

                                        % There was a response, it was in the allowed time window, it was correct

                                        Subject.Negative.TooEarly(nn_Negative) = 0/0;

                                        Subject.Negative.Correct_inCorrect_index(nn_Negative) = 1;
                                        Subject.Negative.NoResponse_index(nn_Negative) = 0/0;

                                        Subject.Negative.Correct_RT(nn_Negative) = round(EEG.event(trial + 1).latency - EEG.event(trial).latency);
                                        Subject.Negative.Incorrect_RT(nn_Negative) = 0/0;

                                    elseif round(EEG.event(trial + 1).latency - EEG.event(trial).latency) < 150

                                        % There was a response, it was too early, it was correct

                                        Subject.Negative.TooEarly(nn_Negative) = 1;

                                        Subject.Negative.Correct_inCorrect_index(nn_Negative) = 0/0;
                                        Subject.Negative.NoResponse_index(nn_Negative) = 0/0;

                                        Subject.Negative.Correct_RT(nn_Negative) = 0/0;
                                        Subject.Negative.Incorrect_RT(nn_Negative) = 0/0;

                                    elseif (round(EEG.event(trial + 1).latency - EEG.event(trial).latency) > 1995)

                                        % There was no response

                                        Subject.Negative.TooEarly(nn_Negative) = 0/0;

                                        Subject.Negative.Correct_inCorrect_index(nn_Negative) = 0/0;
                                        Subject.Negative.NoResponse_index(nn_Negative) = 1;

                                        Subject.Negative.Correct_RT(nn_Negative) = 0/0;
                                        Subject.Negative.Incorrect_RT(nn_Negative) = 0/0;
                                    end

                                elseif strcmp(EEG.event(trial + 1).type, '251') || strcmp(EEG.event(trial + 1).type, 'S251')

                                    if round(EEG.event(trial + 1).latency - EEG.event(trial).latency) >= 150 && round(EEG.event(trial + 1).latency - EEG.event(trial).latency) <= 1995

                                        % There was a response, it was in the allowed time window, it was incorrect

                                        Subject.Negative.TooEarly(nn_Negative) = 0/0;

                                        Subject.Negative.Correct_inCorrect_index(nn_Negative) = 0;
                                        Subject.Negative.NoResponse_index(nn_Negative) = 0/0;

                                        Subject.Negative.Correct_RT(nn_Negative) = round(EEG.event(trial + 1).latency - EEG.event(trial).latency);
                                        Subject.Negative.Incorrect_RT(nn_Negative) = 0/0;

                                    elseif round(EEG.event(trial + 1).latency - EEG.event(trial).latency) < 150

                                        % There was a response, it was too early, it was incorrect

                                        Subject.Negative.TooEarly(nn_Negative) = 1;

                                        Subject.Negative.Correct_inCorrect_index(nn_Negative) = 0/0;
                                        Subject.Negative.NoResponse_index(nn_Negative) = 0/0;

                                        Subject.Negative.Correct_RT(nn_Negative) = 0/0;
                                        Subject.Negative.Incorrect_RT(nn_Negative) = 0/0;

                                    elseif (round(EEG.event(trial + 1).latency - EEG.event(trial).latency) > 1995)

                                        % There was no response

                                        Subject.Negative.TooEarly(nn_Negative) = 0/0;

                                        Subject.Negative.Correct_inCorrect_index(nn_Negative) = 0/0;
                                        Subject.Negative.NoResponse_index(nn_Negative) = 1;

                                        Subject.Negative.Correct_RT(nn_Negative) = 0/0;
                                        Subject.Negative.Incorrect_RT(nn_Negative) = 0/0;

                                    end
                                end

                            elseif strcmp(EEG.event(trial).type, '200') || strcmp(EEG.event(trial).type, '201') || strcmp(EEG.event(trial).type,'S200') || strcmp(EEG.event(trial).type, 'S201')

                                nn_Positive = nn_Positive + 1;

                                if (strcmp(EEG.event(trial + 1).type, '250') || strcmp(EEG.event(trial + 1).type, 'S250'))
                                    if round(EEG.event(trial + 1).latency - EEG.event(trial).latency) >= 150 && round(EEG.event(trial + 1).latency - EEG.event(trial).latency) <= 1995

                                        % There was a response, it was in the allowed time window, it was correct

                                        Subject.Positive.TooEarly(nn_Positive) = 0/0;

                                        Subject.Positive.Correct_inCorrect_index(nn_Positive) = 1;
                                        Subject.Positive.NoResponse_index(nn_Positive) = 0/0;

                                        Subject.Positive.Correct_RT(nn_Positive) = round(EEG.event(trial + 1).latency - EEG.event(trial).latency);
                                        Subject.Positive.Incorrect_RT(nn_Positive) = 0/0;

                                    elseif round(EEG.event(trial + 1).latency - EEG.event(trial).latency) < 150

                                        % There was a response, it was too early, it was correct

                                        Subject.Positive.TooEarly(nn_Positive) = 1;

                                        Subject.Positive.Correct_inCorrect_index(nn_Positive) = 0/0;
                                        Subject.Positive.NoResponse_index(nn_Positive) = 0/0;

                                        Subject.Positive.Correct_RT(nn_Positive) = 0/0;
                                        Subject.Positive.Incorrect_RT(nn_Positive) = 0/0;

                                    elseif (round(EEG.event(trial + 1).latency - EEG.event(trial).latency) > 1995)

                                        % There was no response

                                        Subject.Positive.TooEarly(nn_Positive) = 0/0;

                                        Subject.Positive.Correct_inCorrect_index(nn_Positive) = 0/0;
                                        Subject.Positive.NoResponse_index(nn_Positive) = 1;

                                        Subject.Positive.Correct_RT(nn_Positive) = 0/0;
                                        Subject.Positive.Incorrect_RT(nn_Positive) = 0/0;
                                    end

                                elseif strcmp(EEG.event(trial + 1).type, '251') || strcmp(EEG.event(trial + 1).type, 'S251')

                                    if round(EEG.event(trial + 1).latency - EEG.event(trial).latency) >= 150 && round(EEG.event(trial + 1).latency - EEG.event(trial).latency) <= 1995

                                        % There was a response, it was in the allowed time window, it was incorrect

                                        Subject.Positive.TooEarly(nn_Positive) = 0/0;

                                        Subject.Positive.Correct_inCorrect_index(nn_Positive) = 0;
                                        Subject.Positive.NoResponse_index(nn_Positive) = 0/0;

                                        Subject.Positive.Correct_RT(nn_Positive) = round(EEG.event(trial + 1).latency - EEG.event(trial).latency);
                                        Subject.Positive.Incorrect_RT(nn_Positive) = 0/0;

                                    elseif round(EEG.event(trial + 1).latency - EEG.event(trial).latency) < 150

                                        % There was a response, it was too early, it was incorrect

                                        Subject.Positive.TooEarly(nn_Positive) = 1;

                                        Subject.Positive.Correct_inCorrect_index(nn_Positive) = 0/0;
                                        Subject.Positive.NoResponse_index(nn_Positive) = 0/0;

                                        Subject.Positive.Correct_RT(nn_Positive) = 0/0;
                                        Subject.Positive.Incorrect_RT(nn_Positive) = 0/0;

                                    elseif (round(EEG.event(trial + 1).latency - EEG.event(trial).latency) > 1995)

                                        % There was no response

                                        Subject.Positive.TooEarly(nn_Positive) = 0/0;

                                        Subject.Positive.Correct_inCorrect_index(nn_Positive) = 0/0;
                                        Subject.Positive.NoResponse_index(nn_Positive) = 1;

                                        Subject.Positive.Correct_RT(nn_Positive) = 0/0;
                                        Subject.Positive.Incorrect_RT(nn_Positive) = 0/0;

                                    end
                                end

                            end
                        end

                        save ([inputname(1: end - 4)  '_performance'], 'Subject');

                    catch
                    end
                end
            end
        catch
        end
    end
end

if Generate_performance_table_all_subjects

    % Here, we still pool the entire dataset (including trials that will later be excluded)

    %% Second (group) level analysis

    % First, read the performance data of the subjects that will be analyzed (each group should have no significantly different age)

    for group = 1: length(Group_paper)

        subj = 0;

        try

            cd (['D:\Papers\2025\In preparation\XXX (task microstates)\Data\' Group_paper{group}]);

            a = dir;

            for n = 1: length(a)

                if a(n).isdir && ~strcmp(a(n).name, '.') && ~strcmp(a(n).name, '..')
                    if ((~isempty(find(France_subj_BP_I_Depressed' == str2double(a(n).name(9: end)))) || ~isempty(find(France_subj_BP_I_Euthymic' == str2double(a(n).name(9: end)))) || ~isempty(find(France_subj_BP_II_Depressed' == str2double(a(n).name(9: end)))) || ~isempty(find(France_subj_BP_II_Euthymic' == str2double(a(n).name(9: end)))) || ~isempty(find(France_subj_HC' == str2double(a(n).name(9: end)))) || ~isempty(find(France_subj_Siblings' == str2double(a(n).name(9: end))))))

                        subj = subj + 1;

                        cd (['D:\Papers\2025\In preparation\XXX (task microstates)\Data\' Group_paper{group} '\' a(n).name]);
                        inputname = [a(n).name '_task_processed_performance.mat']; inputpath = pwd;

                        try
                            clear Subject; load(inputname);

                            switch group

                                case 1

                                    France.BP_I.Depressed{subj}.number = Subject.General_info.number;
                                    France.BP_I.Depressed{subj}.age = Subject.General_info.Age_at_examination;
                                    France.BP_I.Depressed{subj}.gender = Subject.General_info.sex;
                                    France.BP_I.Depressed{subj}.med = Subject.meds;
                                    France.BP_I.Depressed{subj}.GAF = Subject.GAF.score;
                                    France.BP_I.Depressed{subj}.MADRS = Subject.MADRS.Apparent_sadness + Subject.MADRS.Reported_sadness + Subject.MADRS.Inner_tension + Subject.MADRS.Reduced_sleep + Subject.MADRS.Reduced_appetite + Subject.MADRS.Concentration_difficulties + Subject.MADRS.Lassitude + Subject.MADRS.Inability_to_feel + Subject.MADRS.Pessimistic_thoughts + Subject.MADRS.Suicidal_thoughts;
                                    France.BP_I.Depressed{subj}.YMRS = Subject.MARS.Elevated_mood + Subject.MARS.Increased_motor_activity_energy + Subject.MARS.Sexual_interest + Subject.MARS.Sleep + Subject.MARS.Irritability + Subject.MARS.Speech + Subject.MARS.Language_thought_disorder + Subject.MARS.Content + Subject.MARS.Disruptive_aggressive_behavior + Subject.MARS.Appearance + Subject.MARS.Insight;

                                    cd ('Neutral');
                                    inputname_data = dir ('*_task_processed_Neutral.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.BP_I.Depressed{subj}.Neutral.TooEarly = Subject.Neutral.TooEarly;
                                    France.BP_I.Depressed{subj}.Neutral.Correct_resp_RT = Subject.Neutral.Correct_RT;
                                    France.BP_I.Depressed{subj}.Neutral.Incorrect_resp_RT = Subject.Neutral.Incorrect_RT;
                                    France.BP_I.Depressed{subj}.Neutral.num_Correct_Incorrect_resp = Subject.Neutral.Correct_inCorrect_index;
                                    France.BP_I.Depressed{subj}.Neutral.num_no_resp = Subject.Neutral.NoResponse_index;
                                    France.BP_I.Depressed{subj}.Neutral.EEG_Data = EEG;

                                    cd ..
                                    cd ('Negative');
                                    inputname_data = dir ('*_task_processed_Negative.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.BP_I.Depressed{subj}.Negative.TooEarly = Subject.Negative.TooEarly;
                                    France.BP_I.Depressed{subj}.Negative.Correct_resp_RT = Subject.Negative.Correct_RT;
                                    France.BP_I.Depressed{subj}.Negative.Incorrect_resp_RT = Subject.Negative.Incorrect_RT;
                                    France.BP_I.Depressed{subj}.Negative.num_Correct_Incorrect_resp = Subject.Negative.Correct_inCorrect_index;
                                    France.BP_I.Depressed{subj}.Negative.num_no_resp = Subject.Negative.NoResponse_index;
                                    France.BP_I.Depressed{subj}.Negative.EEG_Data = EEG;

                                    cd ..
                                    cd ('Positive');
                                    inputname_data = dir ('*_task_processed_Positive.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.BP_I.Depressed{subj}.Positive.TooEarly = Subject.Positive.TooEarly;
                                    France.BP_I.Depressed{subj}.Positive.Correct_resp_RT = Subject.Positive.Correct_RT;
                                    France.BP_I.Depressed{subj}.Positive.Incorrect_resp_RT = Subject.Positive.Incorrect_RT;
                                    France.BP_I.Depressed{subj}.Positive.num_Correct_Incorrect_resp = Subject.Positive.Correct_inCorrect_index;
                                    France.BP_I.Depressed{subj}.Positive.num_no_resp = Subject.Positive.NoResponse_index;
                                    France.BP_I.Depressed{subj}.Positive.EEG_Data = EEG;

                                case 2

                                    France.BP_I.Euthymic{subj}.number =Subject.General_info.number;
                                    France.BP_I.Euthymic{subj}.age = Subject.General_info.Age_at_examination;
                                    France.BP_I.Euthymic{subj}.gender = Subject.General_info.sex;
                                    France.BP_I.Euthymic{subj}.med = Subject.meds;
                                    France.BP_I.Euthymic{subj}.GAF = Subject.GAF.score;
                                    France.BP_I.Euthymic{subj}.MADRS = Subject.MADRS.Apparent_sadness + Subject.MADRS.Reported_sadness + Subject.MADRS.Inner_tension + Subject.MADRS.Reduced_sleep + Subject.MADRS.Reduced_appetite + Subject.MADRS.Concentration_difficulties + Subject.MADRS.Lassitude + Subject.MADRS.Inability_to_feel + Subject.MADRS.Pessimistic_thoughts + Subject.MADRS.Suicidal_thoughts;
                                    France.BP_I.Euthymic{subj}.YMRS = Subject.MARS.Elevated_mood + Subject.MARS.Increased_motor_activity_energy + Subject.MARS.Sexual_interest + Subject.MARS.Sleep + Subject.MARS.Irritability + Subject.MARS.Speech + Subject.MARS.Language_thought_disorder + Subject.MARS.Content + Subject.MARS.Disruptive_aggressive_behavior + Subject.MARS.Appearance + Subject.MARS.Insight;

                                    cd ('Neutral');
                                    inputname_data = dir ('*_task_processed_Neutral.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.BP_I.Euthymic{subj}.Neutral.TooEarly = Subject.Neutral.TooEarly;
                                    France.BP_I.Euthymic{subj}.Neutral.Correct_resp_RT = Subject.Neutral.Correct_RT;
                                    France.BP_I.Euthymic{subj}.Neutral.Incorrect_resp_RT = Subject.Neutral.Incorrect_RT;
                                    France.BP_I.Euthymic{subj}.Neutral.num_Correct_Incorrect_resp = Subject.Neutral.Correct_inCorrect_index;
                                    France.BP_I.Euthymic{subj}.Neutral.num_no_resp = Subject.Neutral.NoResponse_index;
                                    France.BP_I.Euthymic{subj}.Neutral.EEG_Data = EEG;

                                    cd ..
                                    cd ('Negative');
                                    inputname_data = dir ('*_task_processed_Negative.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.BP_I.Euthymic{subj}.Negative.TooEarly = Subject.Negative.TooEarly;
                                    France.BP_I.Euthymic{subj}.Negative.Correct_resp_RT = Subject.Negative.Correct_RT;
                                    France.BP_I.Euthymic{subj}.Negative.Incorrect_resp_RT = Subject.Negative.Incorrect_RT;
                                    France.BP_I.Euthymic{subj}.Negative.num_Correct_Incorrect_resp = Subject.Negative.Correct_inCorrect_index;
                                    France.BP_I.Euthymic{subj}.Negative.num_no_resp = Subject.Negative.NoResponse_index;
                                    France.BP_I.Euthymic{subj}.Negative.EEG_Data = EEG;

                                    cd ..
                                    cd ('Positive');
                                    inputname_data = dir ('*_task_processed_Positive.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.BP_I.Euthymic{subj}.Positive.TooEarly = Subject.Positive.TooEarly;
                                    France.BP_I.Euthymic{subj}.Positive.Correct_resp_RT = Subject.Positive.Correct_RT;
                                    France.BP_I.Euthymic{subj}.Positive.Incorrect_resp_RT = Subject.Positive.Incorrect_RT;
                                    France.BP_I.Euthymic{subj}.Positive.num_Correct_Incorrect_resp = Subject.Positive.Correct_inCorrect_index;
                                    France.BP_I.Euthymic{subj}.Positive.num_no_resp = Subject.Positive.NoResponse_index;
                                    France.BP_I.Euthymic{subj}.Positive.EEG_Data = EEG;

                                case 3

                                    France.BP_II.Depressed{subj}.number =Subject.General_info.number;
                                    France.BP_II.Depressed{subj}.age = Subject.General_info.Age_at_examination;
                                    France.BP_II.Depressed{subj}.gender = Subject.General_info.sex;
                                    France.BP_II.Depressed{subj}.med = Subject.meds;
                                    France.BP_II.Depressed{subj}.GAF = Subject.GAF.score;
                                    France.BP_II.Depressed{subj}.MADRS = Subject.MADRS.Apparent_sadness + Subject.MADRS.Reported_sadness + Subject.MADRS.Inner_tension + Subject.MADRS.Reduced_sleep + Subject.MADRS.Reduced_appetite + Subject.MADRS.Concentration_difficulties + Subject.MADRS.Lassitude + Subject.MADRS.Inability_to_feel + Subject.MADRS.Pessimistic_thoughts + Subject.MADRS.Suicidal_thoughts;
                                    France.BP_II.Depressed{subj}.YMRS = Subject.MARS.Elevated_mood + Subject.MARS.Increased_motor_activity_energy + Subject.MARS.Sexual_interest + Subject.MARS.Sleep + Subject.MARS.Irritability + Subject.MARS.Speech + Subject.MARS.Language_thought_disorder + Subject.MARS.Content + Subject.MARS.Disruptive_aggressive_behavior + Subject.MARS.Appearance + Subject.MARS.Insight;

                                    cd ('Neutral');
                                    inputname_data = dir ('*_task_processed_Neutral.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.BP_II.Depressed{subj}.Neutral.TooEarly = Subject.Neutral.TooEarly;
                                    France.BP_II.Depressed{subj}.Neutral.Correct_resp_RT = Subject.Neutral.Correct_RT;
                                    France.BP_II.Depressed{subj}.Neutral.Incorrect_resp_RT = Subject.Neutral.Incorrect_RT;
                                    France.BP_II.Depressed{subj}.Neutral.num_Correct_Incorrect_resp = Subject.Neutral.Correct_inCorrect_index;
                                    France.BP_II.Depressed{subj}.Neutral.num_no_resp = Subject.Neutral.NoResponse_index;
                                    France.BP_II.Depressed{subj}.Neutral.EEG_Data = EEG;

                                    cd ..
                                    cd ('Negative');
                                    inputname_data = dir ('*_task_processed_Negative.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.BP_II.Depressed{subj}.Negative.TooEarly = Subject.Negative.TooEarly;
                                    France.BP_II.Depressed{subj}.Negative.Correct_resp_RT = Subject.Negative.Correct_RT;
                                    France.BP_II.Depressed{subj}.Negative.Incorrect_resp_RT = Subject.Negative.Incorrect_RT;
                                    France.BP_II.Depressed{subj}.Negative.num_Correct_Incorrect_resp = Subject.Negative.Correct_inCorrect_index;
                                    France.BP_II.Depressed{subj}.Negative.num_no_resp = Subject.Negative.NoResponse_index;
                                    France.BP_II.Depressed{subj}.Negative.EEG_Data = EEG;

                                    cd ..
                                    cd ('Positive');
                                    inputname_data = dir ('*_task_processed_Positive.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.BP_II.Depressed{subj}.Positive.TooEarly = Subject.Positive.TooEarly;
                                    France.BP_II.Depressed{subj}.Positive.Correct_resp_RT = Subject.Positive.Correct_RT;
                                    France.BP_II.Depressed{subj}.Positive.Incorrect_resp_RT = Subject.Positive.Incorrect_RT;
                                    France.BP_II.Depressed{subj}.Positive.num_Correct_Incorrect_resp = Subject.Positive.Correct_inCorrect_index;
                                    France.BP_II.Depressed{subj}.Positive.num_no_resp = Subject.Positive.NoResponse_index;
                                    France.BP_II.Depressed{subj}.Positive.EEG_Data = EEG;

                                case 4

                                    France.BP_II.Euthymic{subj}.number =Subject.General_info.number;
                                    France.BP_II.Euthymic{subj}.age = Subject.General_info.Age_at_examination;
                                    France.BP_II.Euthymic{subj}.gender = Subject.General_info.sex;
                                    France.BP_II.Euthymic{subj}.med = Subject.meds;
                                    France.BP_II.Euthymic{subj}.GAF = Subject.GAF.score;
                                    France.BP_II.Euthymic{subj}.MADRS = Subject.MADRS.Apparent_sadness + Subject.MADRS.Reported_sadness + Subject.MADRS.Inner_tension + Subject.MADRS.Reduced_sleep + Subject.MADRS.Reduced_appetite + Subject.MADRS.Concentration_difficulties + Subject.MADRS.Lassitude + Subject.MADRS.Inability_to_feel + Subject.MADRS.Pessimistic_thoughts + Subject.MADRS.Suicidal_thoughts;
                                    France.BP_II.Euthymic{subj}.YMRS = Subject.MARS.Elevated_mood + Subject.MARS.Increased_motor_activity_energy + Subject.MARS.Sexual_interest + Subject.MARS.Sleep + Subject.MARS.Irritability + Subject.MARS.Speech + Subject.MARS.Language_thought_disorder + Subject.MARS.Content + Subject.MARS.Disruptive_aggressive_behavior + Subject.MARS.Appearance + Subject.MARS.Insight;

                                    cd ('Neutral');
                                    inputname_data = dir ('*_task_processed_Neutral.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.BP_II.Euthymic{subj}.Neutral.TooEarly = Subject.Neutral.TooEarly;
                                    France.BP_II.Euthymic{subj}.Neutral.Correct_resp_RT = Subject.Neutral.Correct_RT;
                                    France.BP_II.Euthymic{subj}.Neutral.Incorrect_resp_RT = Subject.Neutral.Incorrect_RT;
                                    France.BP_II.Euthymic{subj}.Neutral.num_Correct_Incorrect_resp = Subject.Neutral.Correct_inCorrect_index;
                                    France.BP_II.Euthymic{subj}.Neutral.num_no_resp = Subject.Neutral.NoResponse_index;
                                    France.BP_II.Euthymic{subj}.Neutral.EEG_Data = EEG;

                                    cd ..
                                    cd ('Negative');
                                    inputname_data = dir ('*_task_processed_Negative.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.BP_II.Euthymic{subj}.Negative.TooEarly = Subject.Negative.TooEarly;
                                    France.BP_II.Euthymic{subj}.Negative.Correct_resp_RT = Subject.Negative.Correct_RT;
                                    France.BP_II.Euthymic{subj}.Negative.Incorrect_resp_RT = Subject.Negative.Incorrect_RT;
                                    France.BP_II.Euthymic{subj}.Negative.num_Correct_Incorrect_resp = Subject.Negative.Correct_inCorrect_index;
                                    France.BP_II.Euthymic{subj}.Negative.num_no_resp = Subject.Negative.NoResponse_index;
                                    France.BP_II.Euthymic{subj}.Negative.EEG_Data = EEG;

                                    cd ..
                                    cd ('Positive');
                                    inputname_data = dir ('*_task_processed_Positive.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.BP_II.Euthymic{subj}.Positive.TooEarly = Subject.Positive.TooEarly;
                                    France.BP_II.Euthymic{subj}.Positive.Correct_resp_RT = Subject.Positive.Correct_RT;
                                    France.BP_II.Euthymic{subj}.Positive.Incorrect_resp_RT = Subject.Positive.Incorrect_RT;
                                    France.BP_II.Euthymic{subj}.Positive.num_Correct_Incorrect_resp = Subject.Positive.Correct_inCorrect_index;
                                    France.BP_II.Euthymic{subj}.Positive.num_no_resp = Subject.Positive.NoResponse_index;
                                    France.BP_II.Euthymic{subj}.Positive.EEG_Data = EEG;

                                case 5

                                    France.HC{subj}.number =Subject.General_info.number;
                                    France.HC{subj}.age = Subject.General_info.Age_at_examination;
                                    France.HC{subj}.gender = Subject.General_info.sex;
                                    France.HC{subj}.med = Subject.meds;

                                    cd ('Neutral');
                                    inputname_data = dir ('*_task_processed_Neutral.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.HC{subj}.Neutral.TooEarly = Subject.Neutral.TooEarly;
                                    France.HC{subj}.Neutral.Correct_resp_RT = Subject.Neutral.Correct_RT;
                                    France.HC{subj}.Neutral.Incorrect_resp_RT = Subject.Neutral.Incorrect_RT;
                                    France.HC{subj}.Neutral.num_Correct_Incorrect_resp = Subject.Neutral.Correct_inCorrect_index;
                                    France.HC{subj}.Neutral.num_no_resp = Subject.Neutral.NoResponse_index;
                                    France.HC{subj}.Neutral.EEG_Data = EEG;

                                    cd ..
                                    cd ('Negative');
                                    inputname_data = dir ('*_task_processed_Negative.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.HC{subj}.Negative.TooEarly = Subject.Negative.TooEarly;
                                    France.HC{subj}.Negative.Correct_resp_RT = Subject.Negative.Correct_RT;
                                    France.HC{subj}.Negative.Incorrect_resp_RT = Subject.Negative.Incorrect_RT;
                                    France.HC{subj}.Negative.num_Correct_Incorrect_resp = Subject.Negative.Correct_inCorrect_index;
                                    France.HC{subj}.Negative.num_no_resp = Subject.Negative.NoResponse_index;
                                    France.HC{subj}.Negative.EEG_Data = EEG;

                                    cd ..
                                    cd ('Positive');
                                    inputname_data = dir ('*_task_processed_Positive.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.HC{subj}.Positive.TooEarly = Subject.Positive.TooEarly;
                                    France.HC{subj}.Positive.Correct_resp_RT = Subject.Positive.Correct_RT;
                                    France.HC{subj}.Positive.Incorrect_resp_RT = Subject.Positive.Incorrect_RT;
                                    France.HC{subj}.Positive.num_Correct_Incorrect_resp = Subject.Positive.Correct_inCorrect_index;
                                    France.HC{subj}.Positive.num_no_resp = Subject.Positive.NoResponse_index;
                                    France.HC{subj}.Positive.EEG_Data = EEG;

                                case 6

                                    France.Siblings{subj}.number =Subject.General_info.number;
                                    France.Siblings{subj}.age = Subject.General_info.Age_at_examination;
                                    France.Siblings{subj}.gender = Subject.General_info.sex;
                                    France.Siblings{subj}.med = Subject.meds;

                                    cd ('Neutral');
                                    inputname_data = dir ('*_task_processed_Neutral.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.Siblings{subj}.Neutral.TooEarly = Subject.Neutral.TooEarly;
                                    France.Siblings{subj}.Neutral.Correct_resp_RT = Subject.Neutral.Correct_RT;
                                    France.Siblings{subj}.Neutral.Incorrect_resp_RT = Subject.Neutral.Incorrect_RT;
                                    France.Siblings{subj}.Neutral.num_Correct_Incorrect_resp = Subject.Neutral.Correct_inCorrect_index;
                                    France.Siblings{subj}.Neutral.num_no_resp = Subject.Neutral.NoResponse_index;
                                    France.Siblings{subj}.Neutral.EEG_Data = EEG;

                                    cd ..
                                    cd ('Negative');
                                    inputname_data = dir ('*_task_processed_Negative.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.Siblings{subj}.Negative.TooEarly = Subject.Negative.TooEarly;
                                    France.Siblings{subj}.Negative.Correct_resp_RT = Subject.Negative.Correct_RT;
                                    France.Siblings{subj}.Negative.Incorrect_resp_RT = Subject.Negative.Incorrect_RT;
                                    France.Siblings{subj}.Negative.num_Correct_Incorrect_resp = Subject.Negative.Correct_inCorrect_index;
                                    France.Siblings{subj}.Negative.num_no_resp = Subject.Negative.NoResponse_index;
                                    France.Siblings{subj}.Negative.EEG_Data = EEG;

                                    cd ..
                                    cd ('Positive');
                                    inputname_data = dir ('*_task_processed_Positive.mat');
                                    clear EEG; load (inputname_data.name);

                                    France.Siblings{subj}.Positive.TooEarly = Subject.Positive.TooEarly;
                                    France.Siblings{subj}.Positive.Correct_resp_RT = Subject.Positive.Correct_RT;
                                    France.Siblings{subj}.Positive.Incorrect_resp_RT = Subject.Positive.Incorrect_RT;
                                    France.Siblings{subj}.Positive.num_Correct_Incorrect_resp = Subject.Positive.Correct_inCorrect_index;
                                    France.Siblings{subj}.Positive.num_no_resp = Subject.Positive.NoResponse_index;
                                    France.Siblings{subj}.Positive.EEG_Data = EEG;

                            end
                        catch
                            subj = subj - 1;
                            disp('error');
                        end
                    end
                end
            end
        catch

        end
    end

    cd ('D:\Papers\2025\In preparation\XXX (task microstates)\Data\');

    save ('France_data.mat', 'France', '-v7.3');
    clear France

end

if check_Dataset

    % Find how many channels and PCA were omitted

    for group = 1: length(Group_paper)
        try

            cd (['D:\Papers\2025\In preparation\XXX (task microstates)\Data\' Group_paper{group}]);

            a = dir;

            for n = 1: length(a)

                nn_Neutral = 0; nn_Negative = 0; nn_Positive = 0;

                if a(n).isdir && ~strcmp(a(n).name, '.') && ~strcmp(a(n).name, '..')

                    try

                        cd (['D:\Papers\2025\In preparation\XXX (task microstates)\Data\' Group_paper{group} '\' a(n).name]);
                        dd = dir('*CRF.mat');
                        clear Subject; load (dd.name);

                        setFiles = dir(['*processed.set']);
                        setFilenames = {setFiles.name};

                        clear EEG; EEG = pop_loadset('filename', setFilenames, 'filepath', pwd);

                        try
                            num_elec.(Group_paper{group}) = [num_elec.(Group_paper{group}) EEG.nbchan];
                        catch
                            num_elec.(Group_paper{group}) = EEG.nbchan;
                        end

                        try
                            num_ICA.(Group_paper{group}) = [num_ICA.(Group_paper{group}) size(EEG.icaweights, 1)];
                        catch
                            num_ICA.(Group_paper{group}) = size(EEG.icaweights, 1);
                        end

                    end
                end
            end

            Removed_elec.(Group_paper{group}) = 62 - num_elec.(Group_paper{group});
            Removed_ICA_comp.(Group_paper{group}) = num_elec.(Group_paper{group}) - num_ICA.(Group_paper{group});
            Removed_ICA_comp.(Group_paper{group})(find(Removed_ICA_comp.(Group_paper{group}) < 0)) = 0;
        end
    end

    All_removed_elec = [Removed_elec.(Group_paper{1}) Removed_elec.(Group_paper{2}) Removed_elec.(Group_paper{3}) Removed_elec.(Group_paper{4}) Removed_elec.(Group_paper{5}) Removed_elec.(Group_paper{6})]
    All_removed_ICA_comp = [Removed_ICA_comp.(Group_paper{1}) Removed_ICA_comp.(Group_paper{2}) Removed_ICA_comp.(Group_paper{3}) Removed_ICA_comp.(Group_paper{4}) Removed_ICA_comp.(Group_paper{5}) Removed_ICA_comp.(Group_paper{6})]

    % Electrodes

    disp([num2str(100 - (sum(All_removed_elec == 0) / 66) * 100) ' of the datasets required any channel removal']);
    disp(['range across groups ' num2str(min(All_removed_elec)) ' - ' num2str(max(All_removed_elec))]);

    for pop = 1: 6
        disp(['median removed by group: '  num2str(median(Removed_elec.(Group_paper{pop}))) ' in ' Group_paper{pop}]);
    end

    disp(' ');
    disp(' ');

    % ICA components

    disp([num2str(100 - (sum(All_removed_ICA_comp == 0) / 66) * 100) ' of the datasets required any component removal']);
    disp(['range across groups ' num2str(min(All_removed_ICA_comp)) ' - ' num2str(max(All_removed_ICA_comp))]);

    for pop = 1: 6
        disp(['median removed by group: ' num2str(median(Removed_ICA_comp.(Group_paper{pop}))) ' in ' Group_paper{pop}]);
    end
end

if Analyze_data_psychophysics
    load('D:\Papers\2025\In preparation\XXX (task microstates)\Data\France_data.mat');
    Check_demographics(France);
    Check_bad_trials(France);
    Check_RT_normality(France); % Test if RTs are distributed normally
    Check_difference_in_RT_and_accuracy(France); % Non-parametric alternative of 2-way ANOVA (ScheirerRayHare) in case RT are not notmally distributed and Chi square for testing difference in accuracy
end

%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Check_demographics(France)

% Look for normality of age and scores using Shpairo-Wilk test, used for small sample sizes (< 50 samples) and
% check for significant difference across groups using ANOVA1 or its non parametric equivalent - kruskalwallis.

index = 0;

for n = 1: length(France.BP_I.Depressed)
    index = index + 1;
    Age(index) = France.BP_I.Depressed{n}.age;
    Gender(index) = France.BP_I.Depressed{n}.gender;
    GAF(index) =  France.BP_I.Depressed{n}.GAF;
    MADRS(index) =  France.BP_I.Depressed{n}.MADRS;
    YMRS(index) =  France.BP_I.Depressed{n}.YMRS;

    Group_all{index} = 'BP_I_Depressed';

end

index1 = index;

for n = 1: length(France.BP_I.Euthymic)
    index1 = index1 + 1;
    Age(index1) = France.BP_I.Euthymic{n}.age;
    Gender(index1) = France.BP_I.Euthymic{n}.gender;
    GAF(index1) =  France.BP_I.Euthymic{n}.GAF;
    MADRS(index1) =  France.BP_I.Euthymic{n}.MADRS;
    YMRS(index1) =  France.BP_I.Euthymic{n}.YMRS;

    Group_all{index1} = 'BP_I_Euthymic';

end

index2 = index1;

for n = 1: length(France.BP_II.Depressed)
    index2 = index2 + 1;
    Age(index2) = France.BP_II.Depressed{n}.age;
    Gender(index2) = France.BP_II.Depressed{n}.gender;
    GAF(index2) =  France.BP_II.Depressed{n}.GAF;
    MADRS(index2) =  France.BP_II.Depressed{n}.MADRS;
    YMRS(index2) =  France.BP_II.Depressed{n}.YMRS;

    Group_all{index2} = 'BP_II_Depressed';

end

index3 = index2;

for n =1: length(France.BP_II.Euthymic)
    index3 = index3 + 1;
    Age(index3) = France.BP_II.Euthymic{n}.age;
    Gender(index3) = France.BP_II.Euthymic{n}.gender;
    GAF(index3) =  France.BP_II.Euthymic{n}.GAF;
    MADRS(index3) =  France.BP_II.Euthymic{n}.MADRS;
    YMRS(index3) =  France.BP_II.Euthymic{n}.YMRS;

    Group_all{index3} = 'BP_II_Euthymic';

end

index4 = index3;

for n = 1: length(France.HC)
    index4 = index4 + 1;
    Age(index4) = France.HC{n}.age;
    Gender(index4) = France.HC{n}.gender;
    GAF(index4) =  0/0;
    MADRS(index4) =  0/0;
    YMRS(index4) =  0/0;

    Group_all{index4} = 'HC';

end

index5 = index4;

for n = 1: length(France.Siblings)
    index5 = index5 + 1;
    Age(index5) = France.Siblings{n}.age;
    Gender(index5) = France.Siblings{n}.gender;
    GAF(index5) =  0/0;
    MADRS(index5) =  0/0;
    YMRS(index5) =  0/0;

    Group_all{index5} = 'Siblings';

end

disp('------------------------------------------- Demographics ------------------------------------------')
disp(' ');

disp(' _________________Age___________________');
disp(' ');
disp(['BP_I Depressed: ' num2str(round(mean(Age(1: index)) * 10) / 10) ' +- ' num2str(round(std(Age(1: index)) * 10) / 10)]);
disp(' ');
disp('Is age normally distributed?')
disp(' ');
normalitytest(Age(1: index));
disp(' ');
disp(' ');
disp(' ');

disp(['BP_I Euthymic: ' num2str(round(mean(Age(index + 1: index1)) * 10) / 10) ' +- ' num2str(round(std(Age(index + 1: index1)) * 10) / 10)]);
disp(' ');
disp('Is age normally distributed?')
disp(' ');
normalitytest(Age(index + 1: index1));
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Depressed: ' num2str(round(mean(Age(index1 + 1: index2)) * 10) / 10) ' +- ' num2str(round(std(Age(index1 + 1: index2)) * 10) / 10)]);
disp(' ');
disp('Is age normally distributed?')
disp(' ');
normalitytest(Age(index1 + 1: index2));
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Euthymic: '  num2str(round(mean(Age(index2 + 1: index3)) * 10) / 10) ' +- ' num2str(round(std(Age(index2 + 1: index3)) * 10) / 10)]);
disp(' ');
disp('Is age normally distributed?')
disp(' ');
normalitytest(Age(index2 + 1: index3));
disp(' ');
disp(' ');
disp(' ');

disp(['HC: ' num2str(round(mean(Age(index3 + 1: index4)) * 10) / 10) ' +- ' num2str(round(std(Age(index3 + 1: index4)) * 10) / 10)]);
disp(' ');
disp('Is age normally distributed?')
disp(' ');
normalitytest(Age(index3 + 1: index4));
disp(' ');
disp(' ');
disp(' ');

disp(['Siblings: ' num2str(round(mean(Age(index4 + 1: index5)) * 10) / 10) ' +- ' num2str(round(std(Age(index4 + 1: index5)) * 10) / 10)]);
disp(' ');
disp('Is age normally distributed?')
disp(' ');
Results = normalitytest(Age(index4 + 1: index5));
disp(' ');
disp(' ');
disp(' ');

[~, ~, stats] = anova1(Age, Group_all);
[results, ~, ~, gnames] = multcompare(stats);
tbl = array2table(results,"VariableNames", ["Group", "Control Group", "Lower Limit", "Difference", "Upper Limit", "P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

disp(' _________________Gender ___________________');
disp(' ');
disp(['BP_I Depressed: M - ' num2str(length(find(Gender(1: index) == 'M'))) ', F - '  num2str(length(find(Gender(1: index) == 'F')))]);
disp(['BP_I Euthymic: M - ' num2str(length(find(Gender(index + 1: index1) == 'M'))) ', F - '  num2str(length(find(Gender(index + 1: index1) == 'F')))]);
disp(['BP_II Depressed: M - ' num2str(length(find(Gender(index1 + 1: index2) == 'M'))) ', F - '  num2str(length(find(Gender(index1 + 1: index2) == 'F')))]);
disp(['BP_II Euthymic: M - ' num2str(length(find(Gender(index2 + 1: index3) == 'M'))) ', F - '  num2str(length(find(Gender(index2 + 1: index3) == 'F')))]);
disp(['HC: M - ' num2str(length(find(Gender(index3 + 1: index4) == 'M'))) ', F - '  num2str(length(find(Gender(index3 + 1: index4) == 'F')))]);
disp(['Siblings: M - ' num2str(length(find(Gender(index4 + 1: index5) == 'M'))) ', F - '  num2str(length(find(Gender(index4 + 1: index5) == 'F')))]);

disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(' _________________GAF ___________________');
disp(' ');
disp(['BP_I Depressed: ' num2str(round(mean(GAF(1: index)) * 10) / 10) ' +- ' num2str(round(std(GAF(1: index)) * 10) / 10)]);
disp(' ');
disp('Is GAF normally distributed?')
disp(' ');
clear aaa; aaa = GAF(1: index); aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_I Euthymic: ' num2str(round(mean(GAF(index + 1: index1)) * 10) / 10) ' +- ' num2str(round(std(GAF(index + 1: index1)) * 10) / 10)]);
disp(' ');
disp('Is GAF normally distributed?')
disp(' ');
clear aaa; aaa = GAF(index + 1: index1); aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Depressed: ' num2str(round(mean(GAF(index1 + 1: index2)) * 10) / 10) ' +- ' num2str(round(std(GAF(index1 + 1: index2)) * 10) / 10)]);
disp(' ');
disp('Is GAF normally distributed?')
disp(' ');
clear aaa; aaa = GAF(index1 + 1: index2); aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Euthymic: '  num2str(round(mean(GAF(index2 + 1: index3)) * 10) / 10) ' +- ' num2str(round(std(GAF(index2 + 1: index3)) * 10) / 10)]);
disp(' ');
disp('Is GAF normally distributed?')
disp(' ');
clear aaa; aaa = GAF(index2 + 1: index3); aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

clear stats results gnames

[~, ~, stats] = anova1(GAF, Group_all);
[results, ~, ~, gnames] = multcompare(stats);
tbl = array2table(results,"VariableNames", ["Group", "Control Group", "Lower Limit", "Difference", "Upper Limit", "P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(' _________________MADRS ___________________');
disp(' ');
disp(['BP_I Depressed: ' num2str(round(mean(MADRS(1: index)) * 10) / 10) ' +- ' num2str(round(std(MADRS(1: index)) * 10) / 10)]);
disp(' ');
disp('Is MADRS normally distributed?')
disp(' ');
clear aaa; aaa = MADRS(1: index); aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_I Euthymic: ' num2str(round(mean(MADRS(index + 1: index1)) * 10) / 10) ' +- ' num2str(round(std(MADRS(index + 1: index1)) * 10) / 10)]);
disp(' ');
disp('Is MADRS normally distributed?')
disp(' ');
clear aaa; aaa = MADRS(index + 1: index1); aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Depressed: ' num2str(round(mean(MADRS(index1 + 1: index2)) * 10) / 10) ' +- ' num2str(round(std(MADRS(index1 + 1: index2)) * 10) / 10)]);
disp('Is MADRS normally distributed?')
disp(' ');
clear aaa; aaa = MADRS(index1 + 1: index2); aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Euthymic: '  num2str(round(mean(MADRS(index2 + 1: index3)) * 10) / 10) ' +- ' num2str(round(std(MADRS(index2 + 1: index3)) * 10) / 10)]);
disp('Is MADRS normally distributed?')
disp(' ');
clear aaa; aaa = MADRS(index2 + 1: index3); aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

clear stats results gnames

[~, ~, stats] = anova1(MADRS, Group_all);
[results, ~, ~, gnames] = multcompare(stats);
tbl = array2table(results,"VariableNames", ["Group", "Control Group", "Lower Limit", "Difference", "Upper Limit", "P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

disp(' ');
disp(' ');
disp(' ');
disp(' ');
disp(' _________________YMRS ___________________');
disp(' ');
disp(['BP_I Depressed: ' num2str(round(mean(YMRS(1: index)) * 10) / 10) ' +- ' num2str(round(std(YMRS(1: index)) * 10) / 10)]);
disp('Is YMRS normally distributed?')
disp(' ');
clear aaa; aaa = YMRS(1: index); aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_I Euthymic: ' num2str(round(mean(YMRS(index + 1: index1)) * 10) / 10) ' +- ' num2str(round(std(YMRS(index + 1: index1)) * 10) / 10)]);
disp('Is YMRS normally distributed?')
disp(' ');
clear aaa; aaa = YMRS(index + 1: index1); aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Depressed: ' num2str(round(mean(YMRS(index1 + 1: index2)) * 10) / 10) ' +- ' num2str(round(std(YMRS(index1 + 1: index2)) * 10) / 10)]);
disp('Is YMRS normally distributed?')
disp(' ');
clear aaa; aaa = YMRS(index1 + 1: index2); aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

disp(['BP_II Euthymic: '  num2str(round(mean(YMRS(index2 + 1: index3)) * 10) / 10) ' +- ' num2str(round(std(YMRS(index2 + 1: index3)) * 10) / 10)]);
disp('Is YMRS normally distributed?')
disp(' ');
clear aaa; aaa = YMRS(index2 + 1: index3); aaa = aaa(~isnan(aaa));
normalitytest(aaa);
disp(' ');
disp(' ');
disp(' ');

clear stats results gnames

[~, ~, stats] = anova1(YMRS, Group_all);
[results, ~, ~, gnames] = multcompare(stats);
tbl = array2table(results,"VariableNames", ["Group", "Control Group", "Lower Limit", "Difference", "Upper Limit", "P-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Control Group") = gnames(tbl.("Control Group"))

end

% ________________________________________________________________________________________________________________

function Check_bad_trials(France) % Find what is the percentage of no response trials

France_BP_I_Depressed_Neutral_TooEarly = 0; France_BP_I_Depressed_Neutral_NoResponse = 0;
France_BP_I_Depressed_Negative_TooEarly = 0; France_BP_I_Depressed_Negative_NoResponse = 0;
France_BP_I_Depressed_Positive_TooEarly = 0; France_BP_I_Depressed_Positive_NoResponse = 0;

for n = 1: length(France.BP_I.Depressed)
    France_BP_I_Depressed_Neutral_TooEarly = France_BP_I_Depressed_Neutral_TooEarly + length(find(France.BP_I.Depressed{n}.Neutral.TooEarly == 1)) ;
    France_BP_I_Depressed_Neutral_NoResponse = France_BP_I_Depressed_Neutral_NoResponse + length(find(France.BP_I.Depressed{n}.Neutral.num_no_resp == 1)) ;

    France_BP_I_Depressed_Negative_TooEarly = France_BP_I_Depressed_Negative_TooEarly + length(find(France.BP_I.Depressed{n}.Negative.TooEarly == 1)) ;
    France_BP_I_Depressed_Negative_NoResponse = France_BP_I_Depressed_Negative_NoResponse + length(find(France.BP_I.Depressed{n}.Negative.num_no_resp == 1)) ;

    France_BP_I_Depressed_Positive_TooEarly = France_BP_I_Depressed_Positive_TooEarly + length(find(France.BP_I.Depressed{n}.Positive.TooEarly == 1)) ;
    France_BP_I_Depressed_Positive_NoResponse = France_BP_I_Depressed_Positive_NoResponse + length(find(France.BP_I.Depressed{n}.Positive.num_no_resp == 1)) ;
end

France_BP_I_Depressed_Neutral_TooEarly = round(10000 * (France_BP_I_Depressed_Neutral_TooEarly / (40 * length(France.BP_I.Depressed)))) / 100;
France_BP_I_Depressed_Neutral_NoResponse = round(1000 * (France_BP_I_Depressed_Neutral_NoResponse / (40 * length(France.BP_I.Depressed)))) / 100;

France_BP_I_Depressed_Negative_TooEarly = round(10000 * (France_BP_I_Depressed_Negative_TooEarly / (40 * length(France.BP_I.Depressed)))) / 100;
France_BP_I_Depressed_Negative_NoResponse = round(10000 * (France_BP_I_Depressed_Negative_NoResponse / (40 * length(France.BP_I.Depressed)))) / 100;

France_BP_I_Depressed_Positive_TooEarly = round(10000 * (France_BP_I_Depressed_Positive_TooEarly / (40 * length(France.BP_I.Depressed)))) / 100;
France_BP_I_Depressed_Positive_NoResponse = round(10000 * (France_BP_I_Depressed_Positive_NoResponse / (40 * length(France.BP_I.Depressed)))) / 100;

disp('------------------------------------------- Check_bad_trials ------------------------------------------')
disp(' ');
disp(' ------------------------------  BP_I_Depressed --------------------------------');

disp(['BP_I_Depressed_Neutral: Too early: ' num2str(France_BP_I_Depressed_Neutral_TooEarly)]);
disp(['BP_I_Depressed_Neutral: No response: ' num2str(France_BP_I_Depressed_Neutral_NoResponse)]);
disp(' ');

disp(['BP_I_Depressed_Negative: Too early: ' num2str(France_BP_I_Depressed_Negative_TooEarly)]);
disp(['BP_I_Depressed_Negative: No response: ' num2str(France_BP_I_Depressed_Negative_NoResponse)]);
disp(' ');

disp(['BP_I_Depressed_Positive: Too early: ' num2str(France_BP_I_Depressed_Positive_TooEarly)]);
disp(['BP_I_Depressed_Positive: No response: ' num2str(France_BP_I_Depressed_Positive_NoResponse)]);
disp(' ');

France_BP_I_Depressed_Neutral_TooEarly = []; France_BP_I_Depressed_Neutral_NoResponse = [];
France_BP_I_Depressed_Negative_TooEarly = []; France_BP_I_Depressed_Negative_NoResponse = [];
France_BP_I_Depressed_Positive_TooEarly = []; France_BP_I_Depressed_Positive_NoResponse = [];

%%

France_BP_I_Euthymic_Neutral_TooEarly = 0; France_BP_I_Euthymic_Neutral_NoResponse = 0;
France_BP_I_Euthymic_Negative_TooEarly = 0; France_BP_I_Euthymic_Negative_NoResponse = 0;
France_BP_I_Euthymic_Positive_TooEarly = 0; France_BP_I_Euthymic_Positive_NoResponse = 0;

for n = 1: length(France.BP_I.Euthymic)
    France_BP_I_Euthymic_Neutral_TooEarly = France_BP_I_Euthymic_Neutral_TooEarly + length(find(France.BP_I.Euthymic{n}.Neutral.TooEarly == 1)) ;
    France_BP_I_Euthymic_Neutral_NoResponse = France_BP_I_Euthymic_Neutral_NoResponse + length(find(France.BP_I.Euthymic{n}.Neutral.num_no_resp == 1)) ;

    France_BP_I_Euthymic_Negative_TooEarly = France_BP_I_Euthymic_Negative_TooEarly + length(find(France.BP_I.Euthymic{n}.Negative.TooEarly == 1)) ;
    France_BP_I_Euthymic_Negative_NoResponse = France_BP_I_Euthymic_Negative_NoResponse + length(find(France.BP_I.Euthymic{n}.Negative.num_no_resp == 1)) ;

    France_BP_I_Euthymic_Positive_TooEarly = France_BP_I_Euthymic_Positive_TooEarly + length(find(France.BP_I.Euthymic{n}.Positive.TooEarly == 1)) ;
    France_BP_I_Euthymic_Positive_NoResponse = France_BP_I_Euthymic_Positive_NoResponse + length(find(France.BP_I.Euthymic{n}.Positive.num_no_resp == 1)) ;
end

France_BP_I_Euthymic_Neutral_TooEarly = round(10000 * (France_BP_I_Euthymic_Neutral_TooEarly / (40 * length(France.BP_I.Euthymic)))) / 100;
France_BP_I_Euthymic_Neutral_NoResponse = round(1000 * (France_BP_I_Euthymic_Neutral_NoResponse / (40 * length(France.BP_I.Euthymic)))) / 100;

France_BP_I_Euthymic_Negative_TooEarly = round(10000 * (France_BP_I_Euthymic_Negative_TooEarly / (40 * length(France.BP_I.Euthymic)))) / 100;
France_BP_I_Euthymic_Negative_NoResponse = round(10000 * (France_BP_I_Euthymic_Negative_NoResponse / (40 * length(France.BP_I.Euthymic)))) / 100;

France_BP_I_Euthymic_Positive_TooEarly = round(10000 * (France_BP_I_Euthymic_Positive_TooEarly / (40 * length(France.BP_I.Euthymic)))) / 100;
France_BP_I_Euthymic_Positive_NoResponse = round(10000 * (France_BP_I_Euthymic_Positive_NoResponse / (40 * length(France.BP_I.Euthymic)))) / 100;

disp(' ------------------------------  BP_I_Euthymic --------------------------------');

disp(['BP_I_Euthymic_Neutral: Too early: ' num2str(France_BP_I_Euthymic_Neutral_TooEarly)]);
disp(['BP_I_Euthymic_Neutral: No response: ' num2str(France_BP_I_Euthymic_Neutral_NoResponse)]);
disp(' ');

disp(['BP_I_Euthymic_Negative: Too early: ' num2str(France_BP_I_Euthymic_Negative_TooEarly)]);
disp(['BP_I_Euthymic_Negative: No response: ' num2str(France_BP_I_Euthymic_Negative_NoResponse)]);
disp(' ');

disp(['BP_I_Euthymic_Positive: Too early: ' num2str(France_BP_I_Euthymic_Positive_TooEarly)]);
disp(['BP_I_Euthymic_Positive: No response: ' num2str(France_BP_I_Euthymic_Positive_NoResponse)]);
disp(' ');

France.BP_I = [];
France_BP_I_Euthymic_Neutral_TooEarly = []; France_BP_I_Euthymic_Neutral_NoResponse = [];
France_BP_I_Euthymic_Negative_TooEarly = []; France_BP_I_Euthymic_Negative_NoResponse = [];
France_BP_I_Euthymic_Positive_TooEarly = []; France_BP_I_Euthymic_Positive_NoResponse = [];

%%

France_BP_II_Depressed_Neutral_TooEarly = 0; France_BP_II_Depressed_Neutral_NoResponse = 0;
France_BP_II_Depressed_Negative_TooEarly = 0; France_BP_II_Depressed_Negative_NoResponse = 0;
France_BP_II_Depressed_Positive_TooEarly = 0; France_BP_II_Depressed_Positive_NoResponse = 0;

for n = 1: length(France.BP_II.Depressed)
    France_BP_II_Depressed_Neutral_TooEarly = France_BP_II_Depressed_Neutral_TooEarly + length(find(France.BP_II.Depressed{n}.Neutral.TooEarly == 1)) ;
    France_BP_II_Depressed_Neutral_NoResponse = France_BP_II_Depressed_Neutral_NoResponse + length(find(France.BP_II.Depressed{n}.Neutral.num_no_resp == 1)) ;

    France_BP_II_Depressed_Negative_TooEarly = France_BP_II_Depressed_Negative_TooEarly + length(find(France.BP_II.Depressed{n}.Negative.TooEarly == 1)) ;
    France_BP_II_Depressed_Negative_NoResponse = France_BP_II_Depressed_Negative_NoResponse + length(find(France.BP_II.Depressed{n}.Negative.num_no_resp == 1)) ;

    France_BP_II_Depressed_Positive_TooEarly = France_BP_II_Depressed_Positive_TooEarly + length(find(France.BP_II.Depressed{n}.Positive.TooEarly == 1)) ;
    France_BP_II_Depressed_Positive_NoResponse = France_BP_II_Depressed_Positive_NoResponse + length(find(France.BP_II.Depressed{n}.Positive.num_no_resp == 1)) ;
end

France_BP_II_Depressed_Neutral_TooEarly = round(10000 * (France_BP_II_Depressed_Neutral_TooEarly / (40 * length(France.BP_II.Depressed)))) / 100;
France_BP_II_Depressed_Neutral_NoResponse = round(1000 * (France_BP_II_Depressed_Neutral_NoResponse / (40 * length(France.BP_II.Depressed)))) / 100;

France_BP_II_Depressed_Negative_TooEarly = round(10000 * (France_BP_II_Depressed_Negative_TooEarly / (40 * length(France.BP_II.Depressed)))) / 100;
France_BP_II_Depressed_Negative_NoResponse = round(10000 * (France_BP_II_Depressed_Negative_NoResponse / (40 * length(France.BP_II.Depressed)))) / 100;

France_BP_II_Depressed_Positive_TooEarly = round(10000 * (France_BP_II_Depressed_Positive_TooEarly / (40 * length(France.BP_II.Depressed)))) / 100;
France_BP_II_Depressed_Positive_NoResponse = round(10000 * (France_BP_II_Depressed_Positive_NoResponse / (40 * length(France.BP_II.Depressed)))) / 100;

disp(' ------------------------------  BP_II_Depressed --------------------------------');

disp(['BP_II_Depressed_Neutral: Too early: ' num2str(France_BP_II_Depressed_Neutral_TooEarly)]);
disp(['BP_II_Depressed_Neutral: No response: ' num2str(France_BP_II_Depressed_Neutral_NoResponse)]);
disp(' ');

disp(['BP_II_Depressed_Negative: Too early: ' num2str(France_BP_II_Depressed_Negative_TooEarly)]);
disp(['BP_II_Depressed_Negative: No response: ' num2str(France_BP_II_Depressed_Negative_NoResponse)]);
disp(' ');

disp(['BP_II_Depressed_Positive: Too early: ' num2str(France_BP_II_Depressed_Positive_TooEarly)]);
disp(['BP_II_Depressed_Positive: No response: ' num2str(France_BP_II_Depressed_Positive_NoResponse)]);
disp(' ');

France_BP_II_Depressed_Neutral_TooEarly = []; France_BP_II_Depressed_Neutral_NoResponse = [];
France_BP_II_Depressed_Negative_TooEarly = []; France_BP_II_Depressed_Negative_NoResponse = [];
France_BP_II_Depressed_Positive_TooEarly = []; France_BP_II_Depressed_Positive_NoResponse = [];

%%

France_BP_II_Euthymic_Neutral_TooEarly = 0; France_BP_II_Euthymic_Neutral_NoResponse = 0;
France_BP_II_Euthymic_Negative_TooEarly = 0; France_BP_II_Euthymic_Negative_NoResponse = 0;
France_BP_II_Euthymic_Positive_TooEarly = 0; France_BP_II_Euthymic_Positive_NoResponse = 0;

for n = 1: length(France.BP_II.Euthymic)
    France_BP_II_Euthymic_Neutral_TooEarly = France_BP_II_Euthymic_Neutral_TooEarly + length(find(France.BP_II.Euthymic{n}.Neutral.TooEarly == 1)) ;
    France_BP_II_Euthymic_Neutral_NoResponse = France_BP_II_Euthymic_Neutral_NoResponse + length(find(France.BP_II.Euthymic{n}.Neutral.num_no_resp == 1)) ;

    France_BP_II_Euthymic_Negative_TooEarly = France_BP_II_Euthymic_Negative_TooEarly + length(find(France.BP_II.Euthymic{n}.Negative.TooEarly == 1)) ;
    France_BP_II_Euthymic_Negative_NoResponse = France_BP_II_Euthymic_Negative_NoResponse + length(find(France.BP_II.Euthymic{n}.Negative.num_no_resp == 1)) ;

    France_BP_II_Euthymic_Positive_TooEarly = France_BP_II_Euthymic_Positive_TooEarly + length(find(France.BP_II.Euthymic{n}.Positive.TooEarly == 1)) ;
    France_BP_II_Euthymic_Positive_NoResponse = France_BP_II_Euthymic_Positive_NoResponse + length(find(France.BP_II.Euthymic{n}.Positive.num_no_resp == 1)) ;
end

France_BP_II_Euthymic_Neutral_TooEarly = round(10000 * (France_BP_II_Euthymic_Neutral_TooEarly / (40 * length(France.BP_II.Euthymic)))) / 100;
France_BP_II_Euthymic_Neutral_NoResponse = round(1000 * (France_BP_II_Euthymic_Neutral_NoResponse / (40 * length(France.BP_II.Euthymic)))) / 100;

France_BP_II_Euthymic_Negative_TooEarly = round(10000 * (France_BP_II_Euthymic_Negative_TooEarly / (40 * length(France.BP_II.Euthymic)))) / 100;
France_BP_II_Euthymic_Negative_NoResponse = round(10000 * (France_BP_II_Euthymic_Negative_NoResponse / (40 * length(France.BP_II.Euthymic)))) / 100;

France_BP_II_Euthymic_Positive_TooEarly = round(10000 * (France_BP_II_Euthymic_Positive_TooEarly / (40 * length(France.BP_II.Euthymic)))) / 100;
France_BP_II_Euthymic_Positive_NoResponse = round(10000 * (France_BP_II_Euthymic_Positive_NoResponse / (40 * length(France.BP_II.Euthymic)))) / 100;

disp(' ------------------------------  BP_II_Euthymic --------------------------------');

disp(['BP_II_Euthymic_Neutral: Too early: ' num2str(France_BP_II_Euthymic_Neutral_TooEarly)]);
disp(['BP_II_Euthymic_Neutral: No response: ' num2str(France_BP_II_Euthymic_Neutral_NoResponse)]);
disp(' ');

disp(['BP_II_Euthymic_Negative: Too early: ' num2str(France_BP_II_Euthymic_Negative_TooEarly)]);
disp(['BP_II_Euthymic_Negative: No response: ' num2str(France_BP_II_Euthymic_Negative_NoResponse)]);
disp(' ');

disp(['BP_II_Euthymic_Positive: Too early: ' num2str(France_BP_II_Euthymic_Positive_TooEarly)]);
disp(['BP_II_Euthymic_Positive: No response: ' num2str(France_BP_II_Euthymic_Positive_NoResponse)]);
disp(' ');

France.BP_II = [];
France_BP_II_Euthymic_Neutral_TooEarly = []; France_BP_II_Euthymic_Neutral_NoResponse = [];
France_BP_II_Euthymic_Negative_TooEarly = []; France_BP_II_Euthymic_Negative_NoResponse = [];
France_BP_II_Euthymic_Positive_TooEarly = []; France_BP_II_Euthymic_Positive_NoResponse = [];

%%

France_HC_Neutral_TooEarly = 0; France_HC_Neutral_NoResponse = 0;
France_HC_Negative_TooEarly = 0; France_HC_Negative_NoResponse = 0;
France_HC_Positive_TooEarly = 0; France_HC_Positive_NoResponse = 0;

for n = 1: length(France.HC)
    France_HC_Neutral_TooEarly = France_HC_Neutral_TooEarly + length(find(France.HC{n}.Neutral.TooEarly == 1)) ;
    France_HC_Neutral_NoResponse = France_HC_Neutral_NoResponse + length(find(France.HC{n}.Neutral.num_no_resp == 1)) ;

    France_HC_Negative_TooEarly = France_HC_Negative_TooEarly + length(find(France.HC{n}.Negative.TooEarly == 1)) ;
    France_HC_Negative_NoResponse = France_HC_Negative_NoResponse + length(find(France.HC{n}.Negative.num_no_resp == 1)) ;

    France_HC_Positive_TooEarly = France_HC_Positive_TooEarly + length(find(France.HC{n}.Positive.TooEarly == 1)) ;
    France_HC_Positive_NoResponse = France_HC_Positive_NoResponse + length(find(France.HC{n}.Positive.num_no_resp == 1)) ;
end

France_HC_Neutral_TooEarly = round(10000 * (France_HC_Neutral_TooEarly / (40 * length(France.HC)))) / 100;
France_HC_Neutral_NoResponse = round(1000 * (France_HC_Neutral_NoResponse / (40 * length(France.HC)))) / 100;

France_HC_Negative_TooEarly = round(10000 * (France_HC_Negative_TooEarly / (40 * length(France.HC)))) / 100;
France_HC_Negative_NoResponse = round(10000 * (France_HC_Negative_NoResponse / (40 * length(France.HC)))) / 100;

France_HC_Positive_TooEarly = round(10000 * (France_HC_Positive_TooEarly / (40 * length(France.HC)))) / 100;
France_HC_Positive_NoResponse = round(10000 * (France_HC_Positive_NoResponse / (40 * length(France.HC)))) / 100;

disp(' ------------------------------  HC --------------------------------');

disp(['HC_Neutral: Too early: ' num2str(France_HC_Neutral_TooEarly)]);
disp(['HC_Neutral: No response: ' num2str(France_HC_Neutral_NoResponse)]);
disp(' ');

disp(['HC_Negative: Too early: ' num2str(France_HC_Negative_TooEarly)]);
disp(['HC_Negative: No response: ' num2str(France_HC_Negative_NoResponse)]);
disp(' ');

disp(['HC_Positive: Too early: ' num2str(France_HC_Positive_TooEarly)]);
disp(['HC_Positive: No response: ' num2str(France_HC_Positive_NoResponse)]);
disp(' ');

France.HC = [];
France_HC_Neutral_TooEarly = []; France_HC_Neutral_NoResponse = [];
France_HC_Negative_TooEarly = []; France_HC_Negative_NoResponse = [];
France_HC_Positive_TooEarly = []; France_HC_Positive_NoResponse = [];

%%

France_Siblings_Neutral_TooEarly = 0; France_Siblings_Neutral_NoResponse = 0;
France_Siblings_Negative_TooEarly = 0; France_Siblings_Negative_NoResponse = 0;
France_Siblings_Positive_TooEarly = 0; France_Siblings_Positive_NoResponse = 0;

for n = 1: length(France.Siblings)
    France_Siblings_Neutral_TooEarly = France_Siblings_Neutral_TooEarly + length(find(France.Siblings{n}.Neutral.TooEarly == 1)) ;
    France_Siblings_Neutral_NoResponse = France_Siblings_Neutral_NoResponse + length(find(France.Siblings{n}.Neutral.num_no_resp == 1)) ;

    France_Siblings_Negative_TooEarly = France_Siblings_Negative_TooEarly + length(find(France.Siblings{n}.Negative.TooEarly == 1)) ;
    France_Siblings_Negative_NoResponse = France_Siblings_Negative_NoResponse + length(find(France.Siblings{n}.Negative.num_no_resp == 1)) ;

    France_Siblings_Positive_TooEarly = France_Siblings_Positive_TooEarly + length(find(France.Siblings{n}.Positive.TooEarly == 1)) ;
    France_Siblings_Positive_NoResponse = France_Siblings_Positive_NoResponse + length(find(France.Siblings{n}.Positive.num_no_resp == 1)) ;
end

France_Siblings_Neutral_TooEarly = round(10000 * (France_Siblings_Neutral_TooEarly / (40 * length(France.Siblings)))) / 100;
France_Siblings_Neutral_NoResponse = round(1000 * (France_Siblings_Neutral_NoResponse / (40 * length(France.Siblings)))) / 100;

France_Siblings_Negative_TooEarly = round(10000 * (France_Siblings_Negative_TooEarly / (40 * length(France.Siblings)))) / 100;
France_Siblings_Negative_NoResponse = round(10000 * (France_Siblings_Negative_NoResponse / (40 * length(France.Siblings)))) / 100;

France_Siblings_Positive_TooEarly = round(10000 * (France_Siblings_Positive_TooEarly / (40 * length(France.Siblings)))) / 100;
France_Siblings_Positive_NoResponse = round(10000 * (France_Siblings_Positive_NoResponse / (40 * length(France.Siblings)))) / 100;

disp(' ------------------------------  Siblings --------------------------------');

disp(['Siblings_Neutral: Too early: ' num2str(France_Siblings_Neutral_TooEarly)]);
disp(['Siblings_Neutral: No response: ' num2str(France_Siblings_Neutral_NoResponse)]);
disp(' ');

disp(['Siblings_Negative: Too early: ' num2str(France_Siblings_Negative_TooEarly)]);
disp(['Siblings_Negative: No response: ' num2str(France_Siblings_Negative_NoResponse)]);
disp(' ');

disp(['Siblings_Positive: Too early: ' num2str(France_Siblings_Positive_TooEarly)]);
disp(['Siblings_Positive: No response: ' num2str(France_Siblings_Positive_NoResponse)]);
disp(' ');

France.Siblings = [];
France_Siblings_Neutral_TooEarly = []; France_Siblings_Neutral_NoResponse = [];
France_Siblings_Negative_TooEarly = []; France_Siblings_Negative_NoResponse = [];
France_Siblings_Positive_TooEarly = []; France_Siblings_Positive_NoResponse = [];

end

% ____________________________________________________________________________

function Check_RT_normality(France)

% Look for Kolmogorov Test, used for n >= 50

% BP_I Depressed

RT_BP_I_Depressed_Neutral = []; RT_BP_I_Depressed_Negative = []; RT_BP_I_Depressed_Positive = [];
RT_BP_I_Euthymic_Neutral = []; RT_BP_I_Euthymic_Negative = []; RT_BP_I_Euthymic_Positive = [];
RT_BP_II_Depressed_Neutral = []; RT_BP_II_Depressed_Negative = []; RT_BP_II_Depressed_Positive = [];
RT_BP_II_Euthymic_Neutral = []; RT_BP_II_Euthymic_Negative = []; RT_BP_II_Euthymic_Positive = [];
RT_HC_Neutral = []; RT_HC_Negative = []; RT_HC_Positive = [];
RT_Siblings_Neutral = []; RT_Siblings_Negative = []; RT_Siblings_Positive = [];

for n = 1: length(France.BP_I.Depressed)
    RT_BP_I_Depressed_Neutral =[RT_BP_I_Depressed_Neutral France.BP_I.Depressed{n}.Neutral.Correct_resp_RT];
    RT_BP_I_Depressed_Negative = [RT_BP_I_Depressed_Negative France.BP_I.Depressed{n}.Negative.Correct_resp_RT];
    RT_BP_I_Depressed_Positive = [RT_BP_I_Depressed_Positive France.BP_I.Depressed{n}.Positive.Correct_resp_RT];
end

RT_BP_I_Depressed_Neutral = RT_BP_I_Depressed_Neutral(~isnan(RT_BP_I_Depressed_Neutral));
RT_BP_I_Depressed_Negative = RT_BP_I_Depressed_Negative(~isnan(RT_BP_I_Depressed_Negative));
RT_BP_I_Depressed_Positive = RT_BP_I_Depressed_Positive(~isnan(RT_BP_I_Depressed_Positive));

disp(' ');
disp('Are BP_I_Depressed_Neutral RT normally distributed?');
clear H; H = normalitytest(RT_BP_I_Depressed_Neutral);

disp(' ');
disp('Are BP_I_Depressed_Negative RT normally distributed?');
clear H; H = normalitytest(RT_BP_I_Depressed_Negative);

disp(' ');
disp('Are BP_I_Depressed_Positibve RT normally distributed?');
clear H; H = normalitytest(RT_BP_I_Depressed_Positive);

% BP_I Euthymic

for n = 1: length(France.BP_I.Euthymic)
    RT_BP_I_Euthymic_Neutral =[RT_BP_I_Euthymic_Neutral France.BP_I.Euthymic{n}.Neutral.Correct_resp_RT];
    RT_BP_I_Euthymic_Negative = [RT_BP_I_Euthymic_Negative France.BP_I.Euthymic{n}.Negative.Correct_resp_RT];
    RT_BP_I_Euthymic_Positive = [RT_BP_I_Euthymic_Positive France.BP_I.Euthymic{n}.Positive.Correct_resp_RT];
end

disp(' ');
disp(' ');
disp(' ');

RT_BP_I_Euthymic_Neutral = RT_BP_I_Euthymic_Neutral(~isnan(RT_BP_I_Euthymic_Neutral));
RT_BP_I_Euthymic_Negative = RT_BP_I_Euthymic_Negative(~isnan(RT_BP_I_Euthymic_Negative));
RT_BP_I_Euthymic_Positive = RT_BP_I_Euthymic_Positive(~isnan(RT_BP_I_Euthymic_Positive));

disp(' ');
disp('Are BP_I_Euthymic_Neutral RT normally distributed?');
clear H; H = normalitytest(RT_BP_I_Euthymic_Neutral);

disp(' ');
disp('Are BP_I_Euthymic_Negative RT normally distributed?');
clear H; H = normalitytest(RT_BP_I_Euthymic_Negative);

disp(' ');
disp('Are BP_I_Euthymic_Positive RT normally distributed?');
clear H; H = normalitytest(RT_BP_I_Euthymic_Positive);

% BP_II Depressed

for n = 1: length(France.BP_II.Depressed)
    RT_BP_II_Depressed_Neutral =[RT_BP_I_Depressed_Neutral France.BP_II.Depressed{n}.Neutral.Correct_resp_RT];
    RT_BP_II_Depressed_Negative = [RT_BP_I_Depressed_Negative France.BP_II.Depressed{n}.Negative.Correct_resp_RT];
    RT_BP_II_Depressed_Positive = [RT_BP_I_Depressed_Positive France.BP_II.Depressed{n}.Positive.Correct_resp_RT];
end

disp(' ');
disp(' ');
disp(' ');

RT_BP_II_Depressed_Neutral = RT_BP_II_Depressed_Neutral(~isnan(RT_BP_II_Depressed_Neutral));
RT_BP_II_Depressed_Negative = RT_BP_II_Depressed_Negative(~isnan(RT_BP_II_Depressed_Negative));
RT_BP_II_Depressed_Positive = RT_BP_II_Depressed_Positive(~isnan(RT_BP_II_Depressed_Positive));

disp(' ');
disp('Are BP_II_Depressed_Neutral RT normally distributed?');
clear H; H = normalitytest(RT_BP_II_Depressed_Neutral);

disp(' ');
disp('Are BP_II_Depressed_Negative RT normally distributed?');
clear H; H = normalitytest(RT_BP_II_Depressed_Negative);

disp(' ');
disp('Are BP_II_Depressed_Positive RT normally distributed?');
clear H; H = normalitytest(RT_BP_II_Depressed_Positive);

% BP_II Euthymic

for n =1: length(France.BP_II.Euthymic)
    RT_BP_II_Euthymic_Neutral =[RT_BP_I_Euthymic_Neutral France.BP_II.Euthymic{n}.Neutral.Correct_resp_RT];
    RT_BP_II_Euthymic_Negative = [RT_BP_I_Euthymic_Negative France.BP_II.Euthymic{n}.Negative.Correct_resp_RT];
    RT_BP_II_Euthymic_Positive = [RT_BP_I_Euthymic_Positive France.BP_II.Euthymic{n}.Positive.Correct_resp_RT];
end

disp(' ');
disp(' ');
disp(' ');

RT_BP_II_Euthymic_Neutral = RT_BP_II_Euthymic_Neutral(~isnan(RT_BP_II_Euthymic_Neutral));
RT_BP_II_Euthymic_Negative = RT_BP_II_Euthymic_Negative(~isnan(RT_BP_II_Euthymic_Negative));
RT_BP_II_Euthymic_Positive = RT_BP_II_Euthymic_Positive(~isnan(RT_BP_II_Euthymic_Positive));

disp(' ');
disp('Are BP_II_Euthymic_Neutral RT normally distributed?');
clear H; H = normalitytest(RT_BP_II_Euthymic_Neutral);

disp(' ');
disp('Are BP_II_Euthymic_Negative RT normally distributed?');
clear H; H = normalitytest(RT_BP_II_Euthymic_Negative);

disp(' ');
disp('Are BP_II_Euthymic_Positive RT normally distributed?');
clear H; H = normalitytest(RT_BP_II_Euthymic_Positive);

% HC

for n = 1: length(France.HC)
    RT_HC_Neutral =[RT_HC_Neutral France.HC{n}.Neutral.Correct_resp_RT];
    RT_HC_Negative = [RT_HC_Negative France.HC{n}.Negative.Correct_resp_RT];
    RT_HC_Positive = [RT_HC_Positive France.HC{n}.Positive.Correct_resp_RT];
end

disp(' ');
disp(' ');
disp(' ');

RT_HC_Neutral = RT_HC_Neutral(~isnan(RT_HC_Neutral));
RT_HC_Negative = RT_HC_Negative(~isnan(RT_HC_Negative));
RT_HC_Positive = RT_HC_Positive(~isnan(RT_HC_Positive));

disp(' ');
disp('Are HC_Neutral RT normally distributed?');
clear H; H = normalitytest(RT_HC_Neutral);

disp(' ');
disp('Are HC_Negative RT normally distributed?');
clear H; H = normalitytest(RT_HC_Negative);

disp(' ');
disp('Are HC_Positive RT normally distributed?');
clear H; H = normalitytest(RT_HC_Positive);

% Siblings

for n = 1: length(France.Siblings)
    RT_Siblings_Neutral =[RT_Siblings_Neutral France.Siblings{n}.Neutral.Correct_resp_RT];
    RT_Siblings_Negative = [RT_Siblings_Negative France.Siblings{n}.Negative.Correct_resp_RT];
    RT_Siblings_Positive = [RT_Siblings_Positive France.Siblings{n}.Positive.Correct_resp_RT];
end

disp(' ');
disp(' ');
disp(' ');

RT_Siblings_Neutral = RT_Siblings_Neutral(~isnan(RT_Siblings_Neutral));
RT_Siblings_Negative = RT_Siblings_Negative(~isnan(RT_Siblings_Negative));
RT_Siblings_Positive = RT_Siblings_Positive(~isnan(RT_Siblings_Positive));

disp(' ');
disp('Are Siblings_Neutral RT normally distributed?');
clear H; H = normalitytest(RT_Siblings_Neutral);

disp(' ');
disp('Are Siblings_Negative RT normally distributed?');
clear H; H = normalitytest(RT_Siblings_Negative);

disp(' ');
disp('Are Siblings_Positive RT normally distributed?');
clear H; H = normalitytest(RT_Siblings_Positive);

end

function Check_difference_in_RT_and_accuracy(France)

% Incase RT are not normally distributed, use a non-parametric, unbalanced 2-way anova alternative (ScheirerRayHare). If there are differneces,
% use non-parametric 1-way anova (Kruskal-Wallis test) for post-hoc analysis

% Neutral = 100; Negative = 150; Positive = 200;
% 1 - BP_I_Depressed; 2 - BP_I_Eutymic; 3 - BP_II_Depressed; 4 - BP_II_Eutymic; 5 - HC; 6 - Siblings

RT = []; Group = []; Condition = [];

for number_subjects = 1: length(France.BP_I.Depressed)

    RT = [RT France.BP_I.Depressed{number_subjects}.Neutral.Correct_resp_RT];

    for number_trials = 1: length(France.BP_I.Depressed{number_subjects}.Neutral.Correct_resp_RT)
        Group = [Group 1];
        Condition = [Condition 100];
    end

    RT = [RT France.BP_I.Depressed{number_subjects}.Negative.Correct_resp_RT];

    for number_trials = 1: length(France.BP_I.Depressed{number_subjects}.Negative.Correct_resp_RT)
        Group = [Group 1];
        Condition = [Condition 150];
    end

    RT = [RT France.BP_I.Depressed{number_subjects}.Positive.Correct_resp_RT];

    for number_trials = 1: length(France.BP_I.Depressed{number_subjects}.Positive.Correct_resp_RT)
        Group = [Group 1];
        Condition = [Condition 200];
    end
end

for number_subjects = 1: length(France.BP_I.Euthymic)

    RT = [RT France.BP_I.Euthymic{number_subjects}.Neutral.Correct_resp_RT];

    for number_trials = 1: length(France.BP_I.Euthymic{number_subjects}.Neutral.Correct_resp_RT)
        Group = [Group 2];
        Condition = [Condition 100];
    end

    RT = [RT France.BP_I.Euthymic{number_subjects}.Negative.Correct_resp_RT];

    for number_trials = 1: length(France.BP_I.Euthymic{number_subjects}.Negative.Correct_resp_RT)
        Group = [Group 2];
        Condition = [Condition 150];
    end

    RT = [RT France.BP_I.Euthymic{number_subjects}.Positive.Correct_resp_RT];

    for number_trials = 1: length(France.BP_I.Euthymic{number_subjects}.Positive.Correct_resp_RT)
        Group = [Group 2];
        Condition = [Condition 200];
    end
end

for number_subjects = 1: length(France.BP_II.Depressed)

    RT = [RT France.BP_II.Depressed{number_subjects}.Neutral.Correct_resp_RT];

    for number_trials = 1: length(France.BP_II.Depressed{number_subjects}.Neutral.Correct_resp_RT)
        Group = [Group 3];
        Condition = [Condition 100];
    end

    RT = [RT France.BP_II.Depressed{number_subjects}.Negative.Correct_resp_RT];

    for number_trials = 1: length(France.BP_II.Depressed{number_subjects}.Negative.Correct_resp_RT)
        Group = [Group 3];
        Condition = [Condition 150];
    end

    RT = [RT France.BP_II.Depressed{number_subjects}.Positive.Correct_resp_RT];

    for number_trials = 1: length(France.BP_II.Depressed{number_subjects}.Positive.Correct_resp_RT)
        Group = [Group 3];
        Condition = [Condition 200];
    end
end

for number_subjects = 1: length(France.BP_II.Euthymic)

    RT = [RT France.BP_II.Euthymic{number_subjects}.Neutral.Correct_resp_RT];

    for number_trials = 1: length(France.BP_II.Euthymic{number_subjects}.Neutral.Correct_resp_RT)
        Group = [Group 4];
        Condition = [Condition 100];
    end

    RT = [RT France.BP_II.Euthymic{number_subjects}.Negative.Correct_resp_RT];

    for number_trials = 1: length(France.BP_II.Euthymic{number_subjects}.Negative.Correct_resp_RT)
        Group = [Group 4];
        Condition = [Condition 150];
    end

    RT = [RT France.BP_II.Euthymic{number_subjects}.Positive.Correct_resp_RT];

    for number_trials = 1: length(France.BP_II.Euthymic{number_subjects}.Positive.Correct_resp_RT)
        Group = [Group 4];
        Condition = [Condition 200];
    end
end

for number_subjects = 1: length(France.HC)

    RT = [RT France.HC{number_subjects}.Neutral.Correct_resp_RT];

    for number_trials = 1: length(France.HC{number_subjects}.Neutral.Correct_resp_RT)
        Group = [Group 5];
        Condition = [Condition 100];
    end

    RT = [RT France.HC{number_subjects}.Negative.Correct_resp_RT];

    for number_trials = 1: length(France.HC{number_subjects}.Negative.Correct_resp_RT)
        Group = [Group 5];
        Condition = [Condition 150];
    end

    RT = [RT France.HC{number_subjects}.Positive.Correct_resp_RT];

    for number_trials = 1: length(France.HC{number_subjects}.Positive.Correct_resp_RT)
        Group = [Group 5];
        Condition = [Condition 200];
    end
end

for number_subjects = 1: length(France.Siblings)

    RT = [RT France.Siblings{number_subjects}.Neutral.Correct_resp_RT];

    for number_trials = 1: length(France.Siblings{number_subjects}.Neutral.Correct_resp_RT)
        Group = [Group 6];
        Condition = [Condition 100];
    end

    RT = [RT France.Siblings{number_subjects}.Negative.Correct_resp_RT];

    for number_trials = 1: length(France.Siblings{number_subjects}.Negative.Correct_resp_RT)
        Group = [Group 6];
        Condition = [Condition 150];
    end

    RT = [RT France.Siblings{number_subjects}.Positive.Correct_resp_RT];

    for number_trials = 1: length(France.Siblings{number_subjects}.Positive.Correct_resp_RT)
        Group = [Group 6];
        Condition = [Condition 200];
    end
end

Data = [RT', Group', Condition'];

disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
disp(' -                                       Response time                                       -');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
disp(' ');

disp('Effect of group and condition:');
disp(' ');

out = SRH_test(Data, 'Group', 'Condition')

clear RT_p_group RT_t_group RT_stats_group results gnames
[RT_p_group, RT_t_group, RT_stats_group] = kruskalwallis(RT', Group');
[results, ~, ~, gnames] = multcompare(RT_stats_group, 'display', 'off');

disp('Effect of group alone:');
disp(' ');

tbl = array2table(results,"VariableNames", ["Group", "Controling for group", "Lower Limit", "Difference", "Upper Limit", "p-value"]);
tbl.("Group") = gnames(tbl.("Group"));
tbl.("Controling for group") = gnames(tbl.("Controling for group"))

clear RT_p_condition RT_t_condition RT_stats_condition results gnames
[RT_p_condition, RT_t_condition, RT_stats_condition] = kruskalwallis(RT, Condition);
[results, ~, ~, gnames] = multcompare(RT_stats_condition, 'display', 'off');

disp('Effect of condition alone:');
disp(' ');

tbl = array2table(results,"VariableNames", ["Condition", "Controling for condition", "Lower Limit", "Difference", "Upper Limit", "p-value"]);
tbl.("Condition") = gnames(tbl.("Condition"));
tbl.("Controling for condition") = gnames(tbl.("Controling for condition"))

Accuracy = []; Group = []; Condition = [];

% Neutral = 100; Negative = 150; Positive = 200;
% 1 - BP_I_Depressed; 2 - BP_I_Eutymic; 3 - BP_II_Depressed; 4 - BP_II_Eutymic; 5 - HC; 6 - Siblings

for number_subjects = 1: length(France.BP_I.Depressed)

    Accuracy = [Accuracy France.BP_I.Depressed{number_subjects}.Neutral.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.BP_I.Depressed{number_subjects}.Neutral.num_Correct_Incorrect_resp)
        Group = [Group 1];
        Condition = [Condition 100];
    end

    Accuracy = [Accuracy France.BP_I.Depressed{number_subjects}.Negative.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.BP_I.Depressed{number_subjects}.Negative.num_Correct_Incorrect_resp)
        Group = [Group 1];
        Condition = [Condition 150];
    end

    Accuracy = [Accuracy France.BP_I.Depressed{number_subjects}.Positive.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.BP_I.Depressed{number_subjects}.Positive.num_Correct_Incorrect_resp)
        Group = [Group 1];
        Condition = [Condition 200];
    end
end

for number_subjects = 1: length(France.BP_I.Euthymic)

    Accuracy = [Accuracy France.BP_I.Euthymic{number_subjects}.Neutral.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.BP_I.Euthymic{number_subjects}.Neutral.num_Correct_Incorrect_resp)
        Group = [Group 2];
        Condition = [Condition 100];
    end

    Accuracy = [Accuracy France.BP_I.Euthymic{number_subjects}.Negative.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.BP_I.Euthymic{number_subjects}.Negative.num_Correct_Incorrect_resp)
        Group = [Group 2];
        Condition = [Condition 150];
    end

    Accuracy = [Accuracy France.BP_I.Euthymic{number_subjects}.Positive.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.BP_I.Euthymic{number_subjects}.Positive.num_Correct_Incorrect_resp)
        Group = [Group 2];
        Condition = [Condition 200];
    end
end

for number_subjects = 1: length(France.BP_II.Depressed)

    Accuracy = [Accuracy France.BP_II.Depressed{number_subjects}.Neutral.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.BP_II.Depressed{number_subjects}.Neutral.num_Correct_Incorrect_resp)
        Group = [Group 3];
        Condition = [Condition 100];
    end

    Accuracy = [Accuracy France.BP_II.Depressed{number_subjects}.Negative.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.BP_II.Depressed{number_subjects}.Negative.num_Correct_Incorrect_resp)
        Group = [Group 3];
        Condition = [Condition 150];
    end

    Accuracy = [Accuracy France.BP_II.Depressed{number_subjects}.Positive.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.BP_II.Depressed{number_subjects}.Positive.num_Correct_Incorrect_resp)
        Group = [Group 3];
        Condition = [Condition 200];
    end
end

for number_subjects = 1: length(France.BP_II.Euthymic)

    Accuracy = [Accuracy France.BP_II.Euthymic{number_subjects}.Neutral.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.BP_II.Euthymic{number_subjects}.Neutral.num_Correct_Incorrect_resp)
        Group = [Group 4];
        Condition = [Condition 100];
    end

    Accuracy = [Accuracy France.BP_II.Euthymic{number_subjects}.Negative.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.BP_II.Euthymic{number_subjects}.Negative.num_Correct_Incorrect_resp)
        Group = [Group 4];
        Condition = [Condition 150];
    end

    Accuracy = [Accuracy France.BP_II.Euthymic{number_subjects}.Positive.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.BP_II.Euthymic{number_subjects}.Positive.num_Correct_Incorrect_resp)
        Group = [Group 4];
        Condition = [Condition 200];
    end
end

for number_subjects = 1: length(France.HC)

    Accuracy = [Accuracy France.HC{number_subjects}.Neutral.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.HC{number_subjects}.Neutral.num_Correct_Incorrect_resp)
        Group = [Group 5];
        Condition = [Condition 100];
    end

    Accuracy = [Accuracy France.HC{number_subjects}.Negative.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.HC{number_subjects}.Negative.num_Correct_Incorrect_resp)
        Group = [Group 5];
        Condition = [Condition 150];
    end

    Accuracy = [Accuracy France.HC{number_subjects}.Positive.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.HC{number_subjects}.Positive.num_Correct_Incorrect_resp)
        Group = [Group 5];
        Condition = [Condition 200];
    end
end

for number_subjects = 1: length(France.Siblings)

    Accuracy = [Accuracy France.Siblings{number_subjects}.Neutral.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.Siblings{number_subjects}.Neutral.num_Correct_Incorrect_resp)
        Group = [Group 6];
        Condition = [Condition 100];
    end

    Accuracy = [Accuracy France.Siblings{number_subjects}.Negative.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.Siblings{number_subjects}.Negative.num_Correct_Incorrect_resp)
        Group = [Group 6];
        Condition = [Condition 150];
    end

    Accuracy = [Accuracy France.Siblings{number_subjects}.Positive.num_Correct_Incorrect_resp];

    for number_trials = 1: length(France.Siblings{number_subjects}.Positive.num_Correct_Incorrect_resp)
        Group = [Group 6];
        Condition = [Condition 200];
    end
end

clear Accuracy_t_Group Accuracy_chisq_Group Accuracy_p_Group Accuracy_labels_Group Accuracy_t_Condition Accuracy_chisq_Condition Accuracy_p_Condition Accuracy_labels_Condition

[Accuracy_t_Group, Accuracy_chisq_Group, Accuracy_p_Group, Accuracy_labels_Group] = crosstab(Accuracy', Group');
[Accuracy_t_Condition, Accuracy_chisq_Condition, Accuracy_p_Condition, Accuracy_labels_Condition] = crosstab(Accuracy', Condition');

disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
disp(' -                                       Accuracy                                                 -');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
disp(' ');

disp('Effect of group alone:');
disp(' ');
disp(['Chi square = ' num2str(Accuracy_chisq_Group) '. p-value = ' num2str(Accuracy_p_Group)]);

disp(' ');
disp('Effect of condition alone:');
disp(' ');
disp(['Chi square = ' num2str(Accuracy_chisq_Condition) '. p-value = ' num2str(Accuracy_p_Condition)]);

clear GROUP CONDITION Table_Group

GROUP(find(Group == 1)) = "BP\_I\_Depressed";
GROUP(find(Group == 2)) = "BP\_I\_Eutymic";
GROUP(find(Group == 3)) = "BP\_II\_Depressed";
GROUP(find(Group == 4)) = "BP\_II\_Eutymic";
GROUP(find(Group == 5)) = "HC";
GROUP(find(Group == 6)) = "Siblings";

CONDITION(find(Condition == 100)) = "Neutral";
CONDITION(find(Condition == 150)) = "Negative";
CONDITION(find(Condition == 200)) = "Positive";

ACCURACY(find(Accuracy == 0)) = "Wrong";
ACCURACY(find(Accuracy == 1)) = "Correct";

GROUP = GROUP'; CONDITION = CONDITION'; ACCURACY = ACCURACY';
Table_Group = table(GROUP, ACCURACY);Table_Condition = table(CONDITION, ACCURACY);

% For group analysis

Group1_total = sum(Accuracy(Group == 1) ==1) + sum(Accuracy(Group == 1) ==0);
Success_Group1 = sum(Accuracy(Group == 1) ==1);
Percent_Group1 = round((Success_Group1 / Group1_total) * 10000) / 100;

Group2_total = sum(Accuracy(Group == 2) ==1) + sum(Accuracy(Group == 2) ==0);
Success_Group2 = sum(Accuracy(Group == 2) ==1);
Percent_Group2 = round((Success_Group2 / Group2_total) * 10000) / 100;

Group3_total = sum(Accuracy(Group == 3) ==1) + sum(Accuracy(Group == 3) ==0);
Success_Group3 = sum(Accuracy(Group == 3) ==1);
Percent_Group3 = round((Success_Group3 / Group3_total) * 10000) / 100;

Group4_total = sum(Accuracy(Group == 4) ==1) + sum(Accuracy(Group == 4) ==0);
Success_Group4 = sum(Accuracy(Group == 4) ==1);
Percent_Group4 = round((Success_Group4 / Group4_total) * 10000) / 100;

Group5_total = sum(Accuracy(Group == 5) ==1) + sum(Accuracy(Group == 5) ==0);
Success_Group5 = sum(Accuracy(Group == 5) ==1);
Percent_Group5 = round((Success_Group5 / Group5_total) * 10000) / 100;

Group6_total = sum(Accuracy(Group == 6) ==1) + sum(Accuracy(Group == 6) ==0);
Success_Group6 = sum(Accuracy(Group == 6) ==1);
Percent_Group6 = round((Success_Group6 / Group6_total) * 10000) / 100;

% For condition analysis

Cond1_total = sum(Accuracy(Condition == 100) ==1) + sum(Accuracy(Condition == 100) ==0);
Success_Cond1 = sum(Accuracy(Condition == 100) ==1);
Percent_Cond1 = round((Success_Cond1 / Cond1_total) * 10000) / 100;

Cond2_total = sum(Accuracy(Condition == 150) ==1) + sum(Accuracy(Condition == 150) ==0);
Success_Cond2 = sum(Accuracy(Condition == 150) ==1);
Percent_Cond2 = round((Success_Cond2 / Cond2_total) * 10000) / 100;

Cond3_total = sum(Accuracy(Condition == 200) ==1) + sum(Accuracy(Condition == 100) ==0);
Success_Cond3 = sum(Accuracy(Condition == 100) ==1);
Percent_Cond3 = round((Success_Cond1 / Cond3_total) * 10000) / 100;

X_group = [Success_Group1 Group1_total; Success_Group2 Group2_total; Success_Group3 Group3_total; Success_Group4 Group4_total; Success_Group5 Group5_total; Success_Group6 Group6_total];
tmcomptest(X_group, 0.05);

X_condition = [Success_Cond1 Cond1_total; Success_Cond2 Cond2_total; Success_Cond3 Cond3_total];
tmcomptest(X_condition, 0.05);

subplot(1, 2, 1);
bar(1, RT_stats_group.meanranks(1), 'FaceColor',[1 1 1], 'LineWidth', 2); hold on; bar(2.5, RT_stats_group.meanranks(2), 'FaceColor', [1 1 1], 'LineWidth', 2), bar(4, RT_stats_group.meanranks(3), 'FaceColor', [1 1 1], 'LineWidth', 2); bar(5.5, RT_stats_group.meanranks(4), 'FaceColor',[1 1 1], 'LineWidth', 2); bar(7, RT_stats_group.meanranks(5), 'FaceColor',[1 1 1], 'LineWidth', 2); bar(8.5, RT_stats_group.meanranks(6), 'FaceColor',[1 1 1], 'LineWidth', 2);
set(gca, 'XTick', [1 2.5 4 5.5 7 8.5], 'XTickLabel', {'BP\_I\_Depressed', 'BP\_I\_Euthymic', 'BP\_II\_Depressed', 'BP\_II\_Euthymic', 'HC', 'Siblings'}, 'XLim', [0 9], 'LineWidth', 2, 'box', 'off');

subplot(1, 2, 2);
bar(1, RT_stats_condition.meanranks(1), 'FaceColor',[1 1 1], 'LineWidth', 2); hold on; bar(2.5, RT_stats_condition.meanranks(2), 'FaceColor', [1 1 1], 'LineWidth', 2), bar(4, RT_stats_condition.meanranks(3), 'FaceColor', [1 1 1], 'LineWidth', 2);
set(gca, 'XTick', [1 2.5 4], 'XTickLabel', {'Neutral', 'Negative', 'Posiive'}, 'XLim', [0 5], 'LineWidth', 2, 'box', 'off');

set(gcf, 'Color', [1 1 1]);
end
