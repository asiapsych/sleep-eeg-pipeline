function stamp = version_stamp()
% VERSION_STAMP  產生 pipeline 版本戳記(git hash + timestamp + tool versions)
%
% USAGE
%   s = version_stamp();
%   fprintf('Pipeline %s\n', s.label);
%
% 用途:每個 manifest / log / QC PDF 都帶這個戳記,可重現追蹤。
%
% v2.2 ─ 2026/5/27

pipe  = config_pipeline();
paths = config_paths();

stamp = struct();
stamp.pipeline_version = pipe.version;
stamp.timestamp        = datestr(now, 'yyyy-mm-dd HH:MM:SS');
stamp.unix_time        = posixtime(datetime('now'));
stamp.matlab_version   = version();
stamp.host             = char(java.net.InetAddress.getLocalHost.getHostName);

% Git hash(若有 .git 則抓,沒有就 'no_git')
[status, out] = system(sprintf('cd %s && git rev-parse --short HEAD 2>/dev/null', ...
                               paths.PIPELINE_ROOT));
if status == 0
    stamp.git_hash = strtrim(out);
else
    stamp.git_hash = 'no_git';
end

% EEGLAB version
try
    eeglab_v = eeg_getversion();
    stamp.eeglab_version = eeglab_v;
catch
    stamp.eeglab_version = 'unknown';
end

% Plugin 版本(從 plugins 目錄抓)
stamp.plugins = struct();
plugin_dir = fullfile(paths.EEGLAB_DIR, 'plugins');
if exist(plugin_dir, 'dir')
    items = dir(plugin_dir);
    for i = 1:numel(items)
        if items(i).isdir && ~startsWith(items(i).name, '.')
            % 抽取 plugin 名與版本
            nm = items(i).name;
            tok = regexp(nm, '^([A-Za-z_]+)([\d\.]*)$', 'tokens');
            if ~isempty(tok)
                key = tok{1}{1};
                ver = tok{1}{2};
                stamp.plugins.(matlab.lang.makeValidName(key)) = ver;
            else
                stamp.plugins.(matlab.lang.makeValidName(nm)) = '';
            end
        end
    end
end

% Label
stamp.label = sprintf('%s @ %s [git:%s]', ...
    stamp.pipeline_version, stamp.timestamp, stamp.git_hash);
end
