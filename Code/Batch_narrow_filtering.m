% Root directory

root_dir = 'D:\Papers\2026\Submitted\NeuroImage\Data';

% Condition folder names

condition_names = {'Negative', 'Positive', 'Neutral'};

% Start EEGLAB if needed

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Get group folders

group_folders = dir(root_dir);
group_folders = group_folders([group_folders.isdir]);
group_folders = group_folders(~ismember({group_folders.name}, {'.', '..'}));

for g = 1: numel(group_folders)
    group_name = group_folders(g).name;
    group_path = fullfile(root_dir, group_name);
    fprintf('\nProcessing group: %s\n', group_name);

    subject_folders = dir(group_path);
    subject_folders = subject_folders([subject_folders.isdir]);
    subject_folders = subject_folders(~ismember({subject_folders.name}, {'.', '..'}));

    for s = 1: numel(subject_folders)
        subject_name = subject_folders(s).name;
        subject_path = fullfile(group_path, subject_name);
        fprintf('Processing subject: %s\n', subject_name);

        for c = 1: numel(condition_names)
            condition_name = condition_names{c};
            condition_path = fullfile(subject_path, condition_name);

            if ~exist(condition_path, 'dir')
                fprintf('Condition folder not found: %s\n', condition_path);
                continue;
            end

            set_files = dir(fullfile(condition_path, '*.set'));

            if isempty(set_files)
                fprintf('No .set file found in: %s\n', condition_path);
                continue;
            end

            for f = 1: numel(set_files)
                set_name = set_files(f).name;
                [~, base_name, ~] = fileparts(set_name);

                if contains(base_name, '_narrow_filtered')
                    fprintf('Skipping already filtered file: %s\n', set_name);
                    continue;
                end

                new_set_name = [base_name '_narrow_filtered.set'];
                new_mat_name = [base_name '_narrow_filtered.mat'];

                try
                    fprintf('Loading: %s\n', fullfile(condition_path, set_name));
                    EEG = pop_loadset('filename', set_name, 'filepath', condition_path);

                    EEG = pop_eegfiltnew(EEG, 'locutoff', 1, 'hicutoff', 20, 'plotfreqz', 0);
                    EEG = eeg_checkset(EEG);

                    EEG = pop_saveset(EEG, 'filename', new_set_name, 'filepath', condition_path);
                    save(fullfile(condition_path, new_mat_name), 'EEG', '-v7.3');

                    fprintf('Saved: %s\n', fullfile(condition_path, new_set_name));
                    fprintf('Saved: %s\n', fullfile(condition_path, new_mat_name));
                catch ME
                    fprintf('ERROR processing %s\n', fullfile(condition_path, set_name));
                    fprintf('%s\n', ME.message);
                end
            end
        end
    end
end

fprintf('\nAll done.\n');