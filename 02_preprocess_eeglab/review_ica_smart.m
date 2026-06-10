function review_ica_smart(sid_list)
% REVIEW_ICA_SMART  Smart ICA component review(priority-based)
%
% USAGE
%   review_ica_smart();             % 全部 subjects
%   review_ica_smart({'sub01'});
%
% Priority logic(只 review 中信心,自動處理高信心 / 低信心)
%   - 對每 IC 取 max(ICLabel scores) → confidence
%   - confidence >= 0.85: 接受 ICLabel 自動決定(top class)
%   - confidence <= 0.40: 自動 keep(不刪,因為信心太低不確定)
%   - 0.40 < confidence < 0.85: 進 review queue
%
% 操作鍵:
%   r/Enter  → Remove(接受)
%   k        → Keep(覆蓋)
%   s        → Skip subject(全套自動)
%   q        → 儲存進度並結束
%
% 輸出:
%   decisions/<sid>_ica_decisions.json
%     fields: subject, final_remove, auto_remove, reviewer_notes, version_stamp
%
% v2.2 ─ 2026/5/27

if nargin < 1 || isempty(sid_list); sid_list = subjects(); end

paths = config_paths();
thr   = config_thresholds();

fprintf('============================================================\n');
fprintf('  Smart ICA Review\n');
fprintf('  Review window: ICLabel max score in [%.2f, %.2f]\n', ...
        thr.iclabel.review_min, thr.iclabel.review_max);
fprintf('  Outside range = auto-confirm\n');
fprintf('  Keys: r=Remove, k=Keep, s=Skip subj, q=Quit\n');
fprintf('============================================================\n\n');

for i = 1:numel(sid_list)
    sid = sid_list{i};
    ica_path = fullfile(paths.LAB_DATA, [sid '_ica.set']);
    dec_path = fullfile(paths.DECISIONS, [sid '_ica_decisions.json']);

    if ~exist(ica_path, 'file')
        fprintf('[%d/%d] %s: no _ica.set, skip\n', i, numel(sid_list), sid);
        continue;
    end
    if exist(dec_path, 'file')
        fprintf('[%d/%d] %s: decision exists, skip\n', i, numel(sid_list), sid);
        continue;
    end

    fprintf('\n[%d/%d] Loading %s...\n', i, numel(sid_list), sid);
    EEG = pop_loadset('filename', [sid '_ica.set'], 'filepath', paths.LAB_DATA);

    scores  = EEG.etc.ic_classification.ICLabel.classifications;
    classes = EEG.etc.ic_classification.ICLabel.classes;
    n_comp  = size(EEG.icaweights, 1);

    eye_col    = find(strcmpi(classes, 'Eye'));
    muscle_col = find(strcmpi(classes, 'Muscle'));
    brain_col  = find(strcmpi(classes, 'Brain'));

    % 對每 IC 決定 priority
    auto_remove = [];
    auto_keep   = [];
    review_queue = [];
    review_meta = struct('ic', {}, 'top_class', {}, 'top_score', {}, ...
                          'eye', {}, 'muscle', {}, 'brain', {});

    for ci = 1:n_comp
        [top_p, top_idx] = max(scores(ci, :));
        top_class = classes{top_idx};
        eye_s = scores(ci, eye_col);
        mus_s = scores(ci, muscle_col);
        brn_s = scores(ci, brain_col);

        if top_p >= thr.iclabel.review_max
            % 高信心:接受 ICLabel 決定
            if ismember(lower(top_class), {'eye', 'muscle', 'heart', 'line noise', 'channel noise'})
                auto_remove(end+1) = ci; %#ok<AGROW>
            else
                auto_keep(end+1) = ci; %#ok<AGROW>
            end
        elseif top_p <= thr.iclabel.review_min
            % 低信心:預設 keep(不確定就保留)
            auto_keep(end+1) = ci; %#ok<AGROW>
        else
            % 中信心:加入 review queue
            review_queue(end+1) = ci; %#ok<AGROW>
            review_meta(end+1) = struct('ic', ci, 'top_class', top_class, ...
                'top_score', top_p, 'eye', eye_s, 'muscle', mus_s, 'brain', brn_s); %#ok<AGROW>
        end
    end

    fprintf('  Total IC: %d | auto-remove: %d | auto-keep: %d | need review: %d\n', ...
            n_comp, numel(auto_remove), numel(auto_keep), numel(review_queue));

    final_remove   = auto_remove;
    reviewer_notes = {};
    quit_now = false;
    skip_subject = false;

    for q_idx = 1:numel(review_queue)
        ic_num = review_queue(q_idx);
        meta = review_meta(q_idx);

        fprintf('\n  === [%d/%d] IC %d ===\n', q_idx, numel(review_queue), ic_num);
        fprintf('  Top class: %s (%.2f) | Eye %.2f | Muscle %.2f | Brain %.2f\n', ...
                meta.top_class, meta.top_score, meta.eye, meta.muscle, meta.brain);

        % 顯示診斷視窗
        view_ok = false;
        try
            pop_viewprops(EEG, 0, ic_num, [1 70], {}, 0, 'ICLabel');
            drawnow;
            view_ok = true;
        catch err
            fprintf('  pop_viewprops failed (%s), using fallback plot\n', err.message);
        end

        if ~view_ok
            try
                fig = figure('Name', sprintf('%s — IC %d', sid, ic_num), ...
                    'Position', [100 100 1000 700]);
                subplot(2,2,1);
                topoplot(EEG.icawinv(:, ic_num), EEG.chanlocs, ...
                    'electrodes', 'on', 'numcontour', 6);
                title(sprintf('IC %d Scalp', ic_num));
                colorbar;

                subplot(2,2,2);
                act = EEG.icaact(ic_num, :);
                nfft = 2^nextpow2(EEG.srate * 4);
                [pxx, f] = pwelch(act, hann(nfft), nfft/2, nfft, EEG.srate);
                plot(f(f<=70), 10*log10(pxx(f<=70)), 'LineWidth', 1.5);
                xlabel('Frequency (Hz)'); ylabel('Power (dB)');
                title('Power Spectrum');
                xlim([1 70]); grid on;
                xline(4,'--','Delta'); xline(8,'--','Theta');
                xline(13,'--','Alpha'); xline(30,'--','Beta');

                subplot(2,2,[3 4]);
                t_all = (0:length(act)-1) / EEG.srate;
                plot(t_all, act, 'b-', 'LineWidth', 0.5);
                xlabel('Time (s)'); ylabel('Activation');
                title(sprintf('Activation (%d s)', round(t_all(end))));
                grid on;

                sgtitle(sprintf('%s — IC %d | %s %.2f', sid, ic_num, meta.top_class, meta.top_score));
                drawnow;
            catch
            end
        end

        % 等使用者輸入
        while true
            resp = input(sprintf('  [%s IC%d] r=Remove / k=Keep / s=Skip subj / q=Quit: ', ...
                                 sid, ic_num), 's');
            resp = lower(strtrim(resp));
            if isempty(resp); resp = 'r'; end
            if ismember(resp, {'r','k','s','q'}); break; end
            fprintf('  Enter r/k/s/q\n');
        end
        close all;

        switch resp
            case 'q'
                fprintf('  Quitting. Saving partial decision...\n');
                reviewer_notes{end+1} = sprintf('quit_at_IC%d', ic_num); %#ok<AGROW>
                quit_now = true;
                break;
            case 's'
                fprintf('  Skip %s, use auto decisions for rest\n', sid);
                skip_subject = true;
                break;
            case 'r'
                final_remove(end+1) = ic_num; %#ok<AGROW>
                reviewer_notes{end+1} = sprintf('IC%d:remove(reviewed)', ic_num); %#ok<AGROW>
            case 'k'
                reviewer_notes{end+1} = sprintf('IC%d:keep(reviewed)', ic_num); %#ok<AGROW>
        end
    end

    if skip_subject
        % review_queue 中未 review 的 IC 全 keep(預設)
    end

    final_remove = unique(final_remove);

    % 寫 decisions.json
    dec = struct();
    dec.subject        = sid;
    dec.timestamp      = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    dec.final_remove   = final_remove;
    dec.auto_remove    = auto_remove;
    dec.auto_keep      = auto_keep;
    dec.review_queue   = review_queue;
    dec.reviewer_notes = {reviewer_notes};
    dec.version_stamp  = version_stamp();
    fid = fopen(dec_path, 'w');
    fprintf(fid, '%s', jsonencode(dec, 'PrettyPrint', true));
    fclose(fid);

    fprintf('  Saved decision: remove %d IC(s) [%s]\n', ...
            numel(final_remove), num2str(final_remove));
    fprintf('  → %s\n', dec_path);

    if quit_now
        fprintf('\n=== Quit at subject %d/%d. Re-run to continue ===\n', i, numel(sid_list));
        return;
    end
end

close all;
fprintf('\n=== Review done. Next: step2b_apply_ica() ===\n');
end
