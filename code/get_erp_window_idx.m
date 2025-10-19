function idx = get_erp_window_idx(e, S)
% Use ERP.window{e} if present in workspace; otherwise all rows
idx = 1:size(S.MSClass,1);
ws = evalin('base','whos');
hasERP = any(strcmp({ws.name}, 'ERP'));
if hasERP
    ERP = evalin('base','ERP');
    if isstruct(ERP) && isfield(ERP,'window') && numel(ERP.window) >= e && ~isempty(ERP.window{e})
        idx = ERP.window{e};
    end
end
end
