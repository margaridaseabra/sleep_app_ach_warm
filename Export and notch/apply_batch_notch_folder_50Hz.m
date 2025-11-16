function SUMMARY = apply_batch_notch_folder_50Hz(input_dir, varargin)
% Batch 50 Hz notch filtering of EEG stored in .mat files (zero-phase).
%
% For each *.mat:
%   - finds EEG vector (default var 'eeg') and sampling rate (default 'eeg_frequency')
%   - applies 50 Hz notch (and optional harmonics) with filtfilt
%   - saves a new MAT: original variables + eeg_notched + notch_info struct
%
% Usage:
% SUMMARY = notch_folder_50Hz('/path/to/mats', ...
%     'pattern','*.mat', ...
%     'eegVar','eeg', ...
%     'fsVar','eeg_frequency', ...
%     'Q',35, ...
%     'harmonics',true, ...            % notch 50, 100, 150, ... < Nyquist
%     'out_dir','', ...                % default: next to source file
%     'suffix','_notched50.mat', ...
%     'preview_plots',0);              % e.g., 2 => show before/after for first 2 files
%
% Returns: table with per-file status and metadata.

p = inputParser;
addRequired(p,'input_dir',@ischar);
addParameter(p,'pattern','*.mat',@ischar);
addParameter(p,'eegVar','eeg',@ischar);
addParameter(p,'fsVar','eeg_frequency',@ischar);
addParameter(p,'Q',35,@(x)isscalar(x) && x>0);
addParameter(p,'harmonics',true,@islogical);
addParameter(p,'out_dir','',@ischar);
addParameter(p,'suffix','_notched50.mat',@ischar);
addParameter(p,'preview_plots',0,@(x)isscalar(x) && x>=0);
parse(p,input_dir,varargin{:});
S = p.Results;

assert(isfolder(S.input_dir), 'Input folder not found: %s', S.input_dir);
FILES = dir(fullfile(S.input_dir,'**',S.pattern));
FILES = FILES(~[FILES.isdir]);
assert(~isempty(FILES), 'No files matched pattern "%s" under %s', S.pattern, S.input_dir);

SUMMARY = table('Size',[0 9], ...
    'VariableTypes', {'string','double','double','string','string','string','string','logical','string'}, ...
    'VariableNames', {'file','Fs','n_samples','eeg_var','fs_var','out_file','freqs_applied','success','message'});

plot_count = 0;

for i = 1:numel(FILES)
    fpath = fullfile(FILES(i).folder, FILES(i).name);
    [fd,fn,~] = fileparts(fpath);
    try
        D = load(fpath,'-mat');  % load everything (fast enough for metadata)
        % --- find EEG vector
        eeg_field = '';
        if isfield(D, S.eegVar)
            eeg_field = S.eegVar;
        else
            % try some common alternatives
            cand = fieldnames(D);
            hit = find(strcmpi(cand,'eeg') | contains(lower(cand),'eeg'), 1, 'first');
            if ~isempty(hit), eeg_field = cand{hit}; end
        end
        assert(~isempty(eeg_field), 'EEG variable not found. Tried "%s" and similar.', S.eegVar);

        eeg = D.(eeg_field);
        assert(isvector(eeg) && isnumeric(eeg), 'EEG must be a numeric vector.');
        eeg = double(eeg(:));  % force column double
        N = numel(eeg);

        % --- find sampling frequency
        fs_field = '';
        Fs = NaN;
        if isfield(D, S.fsVar)
            fs_field = S.fsVar; Fs = double(D.(fs_field));
        else
            fs_alts = {'fs','Fs','sampling_rate','samplingrate','eeg_fs','eegfrequency','eeg_frequency'};
            for k = 1:numel(fs_alts)
                if isfield(D, fs_alts{k})
                    fs_field = fs_alts{k};
                    Fs = double(D.(fs_field));
                    break;
                end
            end
        end
        assert(isfinite(Fs) && Fs>0, 'Sampling rate not found in expected fields.');

        % --- build notch frequency list
        freqs = 50;
        if S.harmonics
            maxH = floor((Fs/2 - 1) / 50);           % keep a safety margin < Nyquist
            freqs = 50 * (1:maxH);
        end
        freqs = freqs(freqs > 0 & freqs < Fs/2);
        assert(~isempty(freqs), 'No valid notch frequencies below Nyquist for Fs=%.3f', Fs);

        % --- apply sequential iirnotch + filtfilt for each freq
        eeg_f = eeg;
        for f0 = freqs
            w0 = f0/(Fs/2);                  % normalized
            bw = w0/S.Q;                     % bandwidth via Q
            [b,a] = iirnotch(w0, bw);        % requires Signal Processing Toolbox
            eeg_f = filtfilt(b,a, eeg_f);
        end

        % --- save output mat (keep original vars + add eeg_notched & notch_info)
        if isempty(S.out_dir), out_dir = fd; else, out_dir = S.out_dir; end
        if ~isfolder(out_dir), mkdir(out_dir); end
        out_file = fullfile(out_dir, [fn S.suffix]);

        notch_info = struct( ...
            'method','iirnotch+filtfilt', ...
            'Q',S.Q, ...
            'frequencies_Hz',freqs, ...
            'Fs',Fs, ...
            'source_file',fpath, ...
            'eeg_field',eeg_field, ...
            'fs_field',fs_field, ...
            'timestamp',datestr(now) );

        D.eeg_notched = eeg_f;
        D.notch_info  = notch_info;
        save(out_file, '-struct','D','-v7');   % keep v7 (faster) unless very large

        % --- optional quick preview plots for first few files
        if plot_count < S.preview_plots
            plot_count = plot_count + 1;
            L = min(N, min(5*Fs, N));   % plot up to first 5 seconds if long
            t = (0:L-1)/Fs;
            figure('Color','w','Name',['Notch preview — ' fn]);
            subplot(2,1,1); plot(t, eeg(1:L));      xlabel('s'); title(sprintf('%s — raw', fn));
            subplot(2,1,2); plot(t, eeg_f(1:L));    xlabel('s'); title('notched');
        end

        % --- append summary row
        SUMMARY = [SUMMARY; {string(fpath), Fs, N, string(eeg_field), string(fs_field), string(out_file), ...
                             string(strjoin(string(freqs),' ')), true, ""}]; %#ok<AGROW>

    catch ME
        SUMMARY = [SUMMARY; {string(fpath), NaN, NaN, "", "", "", "", false, string(ME.message)}]; %#ok<AGROW>
        warning('notch_folder_50Hz:failed', 'Failed on %s: %s', fpath, ME.message);
        continue;
    end
end

% write a CSV next to input_dir (or out_dir if specified)
summ_dir = ternary(~isempty(S.out_dir), S.out_dir, S.input_dir);
writetable(SUMMARY, fullfile(summ_dir, 'notch_batch_summary.csv'));
fprintf('✅ Done. %d/%d succeeded. Summary: %s\n', nnz(SUMMARY.success), height(SUMMARY), ...
    fullfile(summ_dir, 'notch_batch_summary.csv'));
end

% ----------------- tiny helpers -----------------
function y = ternary(cond, a, b)
y = a; if ~cond, y = b; end
end
