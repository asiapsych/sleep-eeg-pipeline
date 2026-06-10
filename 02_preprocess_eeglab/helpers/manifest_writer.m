function manifest_writer(sid, step_name, info)
% MANIFEST_WRITER  每個 subject 每個 step 的 manifest(JSON)
%
% USAGE
%   info.params = struct('bandpass', [0.5 70], 'bad_ch', {'F7','T7'});
%   info.input  = 'sub01.edf';
%   info.output = 'sub01_step1.set';
%   manifest_writer('sub01', 'step1', info);
%
% 每 subject 累積在 manifests/S<id>.json
%
% v2.2 ─ 2026/5/27

paths = config_paths();
mf_path = fullfile(paths.MANIFEST_DIR, sprintf('%s.json', sid));

% 載入現有 manifest(若有)
if exist(mf_path, 'file')
    fid = fopen(mf_path, 'r');
    raw = fread(fid, inf, '*char')';
    fclose(fid);
    manifest = jsondecode(raw);
else
    manifest = struct();
    manifest.subject = sid;
    manifest.created = datestr(now, 'yyyy-mm-dd HH:MM:SS');
end

% 加入 version stamp(只在第一次)
if ~isfield(manifest, 'version_stamp')
    manifest.version_stamp = version_stamp();
end

% 加入該 step 之資訊
info.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
manifest.(step_name) = info;
manifest.last_updated = datestr(now, 'yyyy-mm-dd HH:MM:SS');

% 寫回 JSON
fid = fopen(mf_path, 'w');
fprintf(fid, '%s', jsonencode(manifest, 'PrettyPrint', true));
fclose(fid);
end
