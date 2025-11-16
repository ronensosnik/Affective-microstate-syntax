function [measure, variant] = parse_filename(f)

% Infer measure and variant from filename

[~, name] = fileparts(f);

if contains(name, 'Coverage', 'IgnoreCase', true)
    measure = 'Coverage';
elseif contains(name, 'Duration', 'IgnoreCase', true)
    measure = 'Duration';
elseif contains(name, 'Occurrence', 'IgnoreCase', true)
    measure = 'Occurrence';
else
    measure = 'Unknown';
end

if contains(name, '_poly', 'IgnoreCase', true)
    variant = 'poly';
else
    variant = 'classes';
end
end