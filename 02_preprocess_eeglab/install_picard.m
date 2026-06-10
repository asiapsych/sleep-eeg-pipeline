% install_picard.m  Auto-install the PICARD plugin for fast ICA.
%
% Requires EEGLAB to already be on the MATLAB path:
%     addpath('/path/to/eeglab'); eeglab nogui;
% Picard is optional — step2_ica falls back to runica if it is unavailable.

eeglab_m = which('eeglab');
if isempty(eeglab_m)
    error(['EEGLAB not found on the MATLAB path. ', ...
           'Run addpath(''/path/to/eeglab''); eeglab nogui; first.']);
end
eeglab_dir = fileparts(eeglab_m);

fprintf('--- Try EEGLAB plugin_askinstall ---\n');
try
    plugin_askinstall('PICARD', [], 1);
    fprintf('plugin_askinstall returned\n');
catch ME
    fprintf('plugin_askinstall failed: %s\n', ME.message);
end

fprintf('\n--- Check picard ---\n');
fprintf('which picard:  %s\n', which('picard'));
fprintf('which Picard:  %s\n', which('Picard'));
fprintf('which runpicard: %s\n', which('runpicard'));

fprintf('\n--- Plugins dir ---\n');
d = dir(fullfile(eeglab_dir, 'plugins'));
for i = 1:numel(d)
    if d(i).isdir && ~startsWith(d(i).name, '.')
        fprintf('  %s\n', d(i).name);
    end
end
