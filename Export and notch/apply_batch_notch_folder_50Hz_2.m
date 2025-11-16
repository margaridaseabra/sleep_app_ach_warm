function SUMMARY = apply_batch_notch_folder_50Hz_2(input_dir, varargin)
% Batch line-noise notch (zero-phase) for EEG inside .mat files.
%
% For each *.mat:
%   - finds EEG vector (default var 'eeg') and sampling rate (default 'eeg_frequency')
%   - applies mains notch at f0 and harmonics (e.g., 50, 100, 150 Hz) using
%     a 2nd-order IIR notch with 3-dB bandwidth specified in **Hz**
%   - saves new MAT with all original vars + eeg_notched + notch_info
%   - writes a summary CSV; optional PSD QA plots
%
% Usage:
% SUMMARY = notch_folder_line_notch('/path/to/mats', ...
%   'pattern','*.mat', ...
%   'eegVar','eeg', ...
%   'fsVar','eeg_frequency', ...
%   'mains',50, ...                % 50 or 60
%   'bw_Hz',3, ...                 % 3 dB bandwidth (Hz) around each harmonic
%   'harmonics',true, ...          % include f0, 2f0, 3f0, ... < Nyquist
%   'out_dir','', ...              % default: next to each source file
%   'suffix','_notched.mat', ...
%   'preview_plots',2, ...         % number of PSD QA PNGs to save
%   'max_plot_freq',100);          % PSD x-limit (Hz)

p = inputParser;
addRequired(p,'input_dir',@ischar);
addParameter(p,'pattern','*.mat',@ischar);
addParameter(p,'eegVar','eeg',@ischar);
addParameter(p,'fsVar','eeg_frequency',@ischar);
addParameter(p,'mains',50,@(x) ismember(x,[50 60]));
addParameter(p,'bw_Hz',3,@(x) isscalar(x) && x>0);
addParameter(p,'harmonics',true,@islogical);
addParameter(p,'out_dir','',@ischar);
addParameter(p,'suffix','_notched.mat',@ischar);
addParameter(p,'preview_plots',0,@(x) isscalar(x) && x>=0);
addParameter(p,'max_plot_freq',100,@(x) isscalar(x) && x>0);
parse(p, input_dir, varargin{:});
S = p.Results;

assert(isfolder(S.input_dir),'Input folder not found: %s',S.input_dir);
FILES = dir(fullfile(S.input_dir,'**',S.pattern));
FILES = FILES(~[FILES.isdir]);
assert(~isempty(FILES),'No files matched "%s" under %s',S.pattern,S.input_dir);

SUMMARY = table('Size',[0 10], ...
    'VariableTypes', {'string','double','double','string','string','string','string','string','logical','string'}, ...
    'VariableNames', {'file','Fs','n_samples','eeg_var','fs_var','out_file','freqs_Hz','bw_Hz','success','message'});

made = 0;

for i = 1:numel(FILES)
    fpath = fullfile(FILES(i).folder, FILES(i).name);
    [fd,fn,~] = fileparts(fpath);

    try
        D = load(fpath,'-mat');

        % --- EEG field
        eeg_field = '';
        if isfield(D,S.eegVar)
            eeg_field = S.eegVar;
        else
            cand = fieldnames(D);
            hit = find(strcmpi(cand,'eeg') | contains(lower(cand),'eeg'),1,'first');
            if ~isempty(hit), eeg_field = cand{hit}; end
        end
        assert(~isempty(eeg_field),'EEG variable not found (looked for "%s").',S.eegVar);
        eeg = D.(eeg_field);  assert(isnumeric(eeg) && isvector(eeg),'EEG must be a numeric vector.');
        eeg = double(eeg(:)); N = numel(eeg);

        % --- Fs field
        fs_field = '';
        Fs = NaN;
        if isfield(D,S.fsVar), fs_field = S.fsVar; Fs = double(D.(fs_field));
        else
            fs_alts = {'fs','Fs','sampling_rate','samplingrate','eeg_fs','eegfrequency','eeg_frequency'};
            for k=1:numel(fs_alts)
                if isfield(D,fs_alts{k}), fs_field = fs_alts{k}; Fs = double(D.(fs_field)); break; end
            end
        end
        assert(isfinite(Fs) && Fs>0,'Sampling rate not found.');

        % --- build f0 list
        if S.harmonics
            maxH = floor((Fs/2 - 1)/S.mains);
            f0_list = S.mains*(1:maxH);
        else
            f0_list = S.mains;
        end
        f0_list = f0_list(f0_list>0 & f0_list<Fs/2);
        assert(~isempty(f0_list),'No valid notch freqs below Nyquist.');

        % --- apply notch (Hz bandwidth)
        eeg_f = notch_multi(eeg, Fs, f0_list, S.bw_Hz);

        % --- save output
        out_dir = fd; if ~isempty(S.out_dir), out_dir = S.out_dir; end
        if ~isfolder(out_dir), mkdir(out_dir); end
        out_file = fullfile(out_dir, [fn S.suffix]);

        notch_info = struct( ...
            'method','iirnotch+filtfilt', ...
            'mains',S.mains, ...
            'frequencies_Hz',f0_list, ...
            'bw_Hz',S.bw_Hz, ...
            'Fs',Fs, ...
            'eeg_field',eeg_field, ...
            'fs_field',fs_field, ...
            'source_file',fpath, ...
            'timestamp',datestr(now));
        D.eeg_notched = eeg_f;
        D.notch_info  = notch_info;
        save(out_file,'-struct','D','-v7');

        % --- QA PSD plot (optional)
        if made < S.preview_plots
            made = made + 1;
            qa_psd_plot(eeg, eeg_f, Fs, f0_list, S.max_plot_freq, fullfile(out_dir, [fn '_notch_QA.png']));
        end

        % --- summary
        SUMMARY = [SUMMARY; {string(fpath), Fs, N, string(eeg_field), string(fs_field), ...
                             string(out_file), string(strjoin(string(f0_list),' ')), ...
                             string(num2str(S.bw_Hz)), true, ""}]; %#ok<AGROW>

    catch ME
        SUMMARY = [SUMMARY; {string(fpath), NaN, NaN, "", "", "", "", string(num2str(S.bw_Hz)), false, string(ME.message)}]; %#ok<AGROW>
        warning('notch:failed','Failed on %s: %s', fpath, ME.message);
    end
end

% write CSV summary
summ_dir = S.out_dir; if isempty(summ_dir), summ_dir = S.input_dir; end
writetable(SUMMARY, fullfile(summ_dir, sprintf('notch_batch_summary_%dHz_bw%.1f.csv', S.mains, S.bw_Hz)));
fprintf('✅ Notch done: %d/%d OK. CSV: %s\n', nnz(SUMMARY.success), height(SUMMARY), ...
    fullfile(summ_dir, sprintf('notch_batch_summary_%dHz_bw%.1f.csv', S.mains, S.bw_Hz)));
end

% ---------- helpers ----------
function x_out = notch_multi(x_in, Fs, f0_list, bw_Hz)
x = double(x_in(:));
x = x - mean(x,'omitnan');          % de-mean to reduce DC leakage
for f0 = f0_list(:).'
    if f0 >= Fs/2 || f0<=0, continue; end
    w0 = f0/(Fs/2);                  % normalized center
    bw = bw_Hz/(Fs/2);               % normalized 3-dB bandwidth
    [b,a] = iirnotch(w0, bw);        % 2nd-order notch
    x = filtfilt(b,a,x);             % zero-phase
end
x_out = reshape(x, size(x_in));
end

function qa_psd_plot(raw, notched, Fs, f0_list, fmax, out_png)
nw = round(8*Fs);                    % 4 s Welch window
nover = round(0.5*nw);
nfft  = 2^nextpow2(ceil(Fs/0.25));

fgrid = 0:0.5:fmax;
% [Pr,fr] = pwelch(raw, nw, nover, fgrid, Fs, 'onesided');
% [Pn,fn] = pwelch(notched, nw, nover, fgrid, Fs, 'onesided');
[Pr,fr] = pwelch(raw, nw, nover, nfft, Fs, 'onesided');
[Pn,fn] = pwelch(notched, nw, nover, nfft, Fs, 'onesided');

fig = figure('Color','w','Visible','off','Position',[100 100 1000 500]);
plot(fr,Pr,'Color',[0.75 0.75 0.9]); hold on
plot(fn,Pn,'b','LineWidth',1.2);
for f0 = f0_list
    xline(f0,'r--','LineWidth',1);
end
xlim([0 fmax]); ylabel('PSD'); xlabel('Frequency (Hz)');
legend({'raw','notched'},'Location','northeast'); box off
title(sprintf('Line notch QA — f0=[%s] Hz, bw=%.1f Hz', strjoin(string(f0_list),' '), f0_list*0+NaN)); %#ok<NASGU>
% re-do title without NaN hack:
title(sprintf('Line notch QA — f0=[%s] Hz, bw=%.1f Hz', strjoin(string(f0_list),' '), f0_list(1)*0+0 + 0 + 0 + 0 + 0 + 0 + 0 + 0)); % dummy math to keep old MATLAB happy
title(sprintf('Line notch QA — f0=[%s] Hz', strjoin(string(f0_list),' ')));
exportgraphics(fig, out_png, 'Resolution',150);
close(fig);
end
