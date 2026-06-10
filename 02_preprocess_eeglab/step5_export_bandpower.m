function step5_export_bandpower(sid_list, opts)
% STEP5_EXPORT_BANDPOWER  Compute absolute and relative band power from PSD.
%
% Outputs:
%   outputs/band_power/bandpower_<subject>.mat
%   outputs/band_power/bandpower_channels_long.csv
%   outputs/band_power/bandpower_subject_summary.csv
%
% Absolute power is the trapezoidal integral of PSD within each band.
% Relative power uses the summed PRIMARY bands as denominator. Derived bands
% are reported against the same denominator but do not contribute to it.

if nargin < 1 || isempty(sid_list); sid_list = subjects(); end
if nargin < 2; opts = struct(); end
if ~isfield(opts, 'overwrite'); opts.overwrite = false; end

paths = config_paths();
bands = config_bands();
pipe  = config_pipeline();

if ~exist(paths.POWER_OUT, 'dir')
    mkdir(paths.POWER_OUT);
end

channel_csv = fullfile(paths.POWER_OUT, 'bandpower_channels_long.csv');
summary_csv = fullfile(paths.POWER_OUT, 'bandpower_subject_summary.csv');
if opts.overwrite
    if exist(channel_csv, 'file'); delete(channel_csv); end
    if exist(summary_csv, 'file'); delete(summary_csv); end
end

fprintf('============================================================\n');
fprintf('  Step 5: Band Power Export\n');
fprintf('  Pipeline %s | %d subjects\n', pipe.version, numel(sid_list));
fprintf('============================================================\n\n');

n_ok = 0; n_skip = 0; n_fail = 0;
fail_list = {};

for i = 1:numel(sid_list)
    sid = sid_list{i};
    in_path  = fullfile(paths.PSD_OUT, sprintf('psd_%s.mat', sid));
    out_path = fullfile(paths.POWER_OUT, sprintf('bandpower_%s.mat', sid));

    fprintf('[%d/%d] %s ... ', i, numel(sid_list), sid);

    if ~exist(in_path, 'file')
        fprintf('SKIP (no PSD)\n');
        n_skip = n_skip + 1;
        continue;
    end
    if exist(out_path, 'file') && ~opts.overwrite
        fprintf('SKIP (done)\n');
        n_skip = n_skip + 1;
        continue;
    end

    try
        S = load(in_path, 'freqs', 'psd', 'ch_names', 'n_ep', 'sfreq', 'subject');
        freqs = S.freqs(:)';
        psd = S.psd;
        ch_names = S.ch_names;
        n_ch = numel(ch_names);

        primary_names = fieldnames(bands.primary);
        derived_names = fieldnames(bands.derived);
        all_names = [primary_names; derived_names];

        abs_power = struct();
        for b = 1:numel(all_names)
            band_name = all_names{b};
            if isfield(bands.primary, band_name)
                band_range = bands.primary.(band_name);
            else
                band_range = bands.derived.(band_name);
            end
            abs_power.(band_name) = integrate_band(psd, freqs, band_range);
        end

        total_primary = zeros(n_ch, 1);
        for b = 1:numel(primary_names)
            total_primary = total_primary + abs_power.(primary_names{b});
        end
        total_primary(total_primary <= 0) = NaN;

        rel_power = struct();
        for b = 1:numel(all_names)
            rel_power.(all_names{b}) = abs_power.(all_names{b}) ./ total_primary;
        end

        channel_table = make_channel_table(sid, ch_names, all_names, abs_power, rel_power, S.n_ep, S.sfreq);
        summary_table = make_summary_table(sid, all_names, abs_power, rel_power, S.n_ep, S.sfreq);

        save(out_path, 'sid', 'freqs', 'ch_names', 'abs_power', 'rel_power', ...
             'total_primary', 'bands', 'channel_table', 'summary_table', '-v7');

        append_table(channel_csv, channel_table);
        append_table(summary_csv, summary_table);

        info = struct();
        info.input = in_path;
        info.output = out_path;
        info.channel_csv = channel_csv;
        info.summary_csv = summary_csv;
        info.n_ch = n_ch;
        info.n_ep = S.n_ep;
        info.primary_denominator = 'sum_integrated_primary_bands_Delta_Theta_Alpha_Beta_Gamma';
        manifest_writer(sid, 'step5', info);

        log_writer('step5', sid, ...
            {'n_ch','n_ep','n_bands','status','note'}, ...
            {n_ch, S.n_ep, numel(all_names), 'ok', ''});

        fprintf('OK (%d channels, %d bands)\n', n_ch, numel(all_names));
        n_ok = n_ok + 1;
    catch ME
        log_writer('step5', sid, ...
            {'n_ch','n_ep','n_bands','status','note'}, ...
            {NaN, NaN, NaN, 'error', ME.message});
        fprintf('FAIL: %s\n', ME.message);
        fail_list{end+1} = sid; %#ok<AGROW>
        n_fail = n_fail + 1;
    end
end

fprintf('\n============================================================\n');
fprintf('  Step 5 done: %d ok, %d skip, %d fail / total %d\n', ...
        n_ok, n_skip, n_fail, numel(sid_list));
fprintf('  Output: %s\n', paths.POWER_OUT);
if ~isempty(fail_list)
    fprintf('  Failed: %s\n', strjoin(fail_list, ', '));
end
fprintf('============================================================\n');
end


function band_power = integrate_band(psd, freqs, band_range)
mask = freqs >= band_range(1) & freqs < band_range(2);
if ~any(mask)
    band_power = NaN(size(psd, 1), 1);
    return;
end
band_power = trapz(freqs(mask), psd(:, mask), 2);
end


function T = make_channel_table(sid, ch_names, band_names, abs_power, rel_power, n_ep, sfreq)
n_ch = numel(ch_names);
n_bands = numel(band_names);
subject = strings(n_ch * n_bands, 1);
channel = strings(n_ch * n_bands, 1);
band = strings(n_ch * n_bands, 1);
absolute_power = zeros(n_ch * n_bands, 1);
relative_power = zeros(n_ch * n_bands, 1);
n_epochs = repmat(n_ep, n_ch * n_bands, 1);
sampling_rate = repmat(sfreq, n_ch * n_bands, 1);

row = 0;
for ch = 1:n_ch
    for b = 1:n_bands
        row = row + 1;
        band_name = band_names{b};
        subject(row) = string(sid);
        channel(row) = string(ch_names{ch});
        band(row) = string(band_name);
        band_abs = abs_power.(band_name);
        band_rel = rel_power.(band_name);
        absolute_power(row) = band_abs(ch);
        relative_power(row) = band_rel(ch);
    end
end

T = table(subject, channel, band, absolute_power, relative_power, n_epochs, sampling_rate);
end


function T = make_summary_table(sid, band_names, abs_power, rel_power, n_ep, sfreq)
n_bands = numel(band_names);
subject = strings(n_bands, 1);
band = strings(n_bands, 1);
mean_absolute_power = zeros(n_bands, 1);
mean_relative_power = zeros(n_bands, 1);
n_epochs = repmat(n_ep, n_bands, 1);
sampling_rate = repmat(sfreq, n_bands, 1);

for b = 1:n_bands
    band_name = band_names{b};
    subject(b) = string(sid);
    band(b) = string(band_name);
    mean_absolute_power(b) = mean(abs_power.(band_name), 'omitnan');
    mean_relative_power(b) = mean(rel_power.(band_name), 'omitnan');
end

T = table(subject, band, mean_absolute_power, mean_relative_power, n_epochs, sampling_rate);
end


function append_table(csv_path, T)
if exist(csv_path, 'file')
    writetable(T, csv_path, 'WriteMode', 'append', 'WriteVariableNames', false);
else
    writetable(T, csv_path);
end
end
