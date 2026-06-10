function log_writer(step_name, sid, fields, values)
% LOG_WRITER  統一 log 寫入(每 step 一個 CSV)
%
% USAGE
%   log_writer('step1', 'sub01', ...
%       {'n_chan_in','srate','duration_s','bad_channels','n_bad','status','note'}, ...
%       {19, 300, 350.5, 'F7;T7', 2, 'ok', ''})
%
% 第一次呼叫某 step 會建 csv + header;之後 append。
%
% v2.2 ─ 2026/5/27

paths = config_paths();
log_path = fullfile(paths.LOGS, sprintf('%s_log.csv', step_name));

% 若檔不存在,寫 header
if ~exist(log_path, 'file')
    fid = fopen(log_path, 'w');
    fprintf(fid, 'subject,%s\n', strjoin(fields, ','));
    fclose(fid);
end

% Append 該筆
fid = fopen(log_path, 'a');
fprintf(fid, '%s', sid);
for i = 1:numel(values)
    v = values{i};
    if islogical(v)
        fprintf(fid, ',%d', double(v));
    elseif isnumeric(v)
        if isempty(v)
            fprintf(fid, ',');
        elseif isnan(v)
            fprintf(fid, ',');
        elseif mod(v, 1) == 0
            fprintf(fid, ',%d', v);
        else
            fprintf(fid, ',%.3f', v);
        end
    else
        % 字串(去除逗號避免破壞 csv)
        s = strrep(char(v), ',', ';');
        s = strrep(s, '"', '''');
        if any(s == newline) || any(s == sprintf('\n')) || any(s == sprintf('\r'))
            s = strrep(s, newline, ' ');
        end
        fprintf(fid, ',"%s"', s);
    end
end
fprintf(fid, '\n');
fclose(fid);
end
