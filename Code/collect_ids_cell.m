    function ids = collect_ids_cell(cellarr)
        ids = NaN(1, numel(cellarr));
        for k = 1: numel(cellarr)
            s = cellarr{k};
            if isempty(s) || ~isfield(s, 'number')
                ids(k) = NaN;
            else
                ids(k) = s.number;
            end
        end
    end
