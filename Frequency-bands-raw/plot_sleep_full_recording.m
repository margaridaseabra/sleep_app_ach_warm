function plot_sleep_full_recording()
% plot_sleep_full_recording.m
% MATLAB version of your Python script:
% - Load .mat signals (eeg/emg/eeg_frequency)
% - Load 1 Hz scoring CSV
% - Load warming intervals from XLSX
% - Compute delta/sigma band power via Hilbert method (1 Hz bins)
% - Plot EEG + delta + sigma + warming strip
% - Batch export to multi-page PDF (N per group)
%
% Requires: Signal Processing Toolbox (butter, filtfilt, hilbert)

% ===========================
% CONFIG
% ===========================
SIGNAL_DIR  = "/Users/margaridaseabra/24.11 signalnotscored";
SCORES_DIR  = "/Users/margaridaseabra/24.11scores";
WARMING_XLSX = "/Users/margaridaseabra/sleep_app_ach_warm/Ambtemp prep-3h-selection/ambtemp-times-all-mice.xlsx";

DELTA_BAND = [0.5, 5.0];
SIGMA_BAND = [10.0, 15.0];

POWER_BIN_SEC    = 1.0;  % align with 1 Hz scoring
POWER_SMOOTH_SEC = 5.0;  % set 0 to disable
EEG_MAX_PLOT_POINTS = 250000;

% State colors (RGBA-ish: MATLAB uses FaceAlpha separately)
STATE_COLORS = containers.Map('KeyType','double','ValueType','any');
STATE_COLORS(0)  = [0.90 0.80 0.95 0.35]; % Wake
STATE_COLORS(1)  = [1.00 0.60 0.60 0.35]; % NREM
STATE_COLORS(2)  = [0.60 1.00 0.60 0.35]; % REM
STATE_COLORS(15) = [0.70 0.70 1.00 0.35]; % MA

STATE_NAMES = containers.Map('KeyType','double','ValueType','char');
STATE_NAMES(0)  = "Wake";
STATE_NAMES(1)  = "NREM";
STATE_NAMES(2)  = "REM";
STATE_NAMES(15) = "MA";

% crude genotype inference from filename (edit to match your naming)
GROUP_PATTERNS = struct();
GROUP_PATTERNS.WT  = "(^|[\W_])WT($|[\W_])|_WT_|WT_";
GROUP_PATTERNS.APP = "(^|[\W_])APP($|[\W_])|_APP_|APP_";

% ===========================
% MAIN MENU
% ===========================
mats = dir(fullfile(SIGNAL_DIR, "*.mat"));
if isempty(mats)
    fprintf("No .mat files found in %s\n", SIGNAL_DIR);
    return;
end

fprintf("Found %d .mat files.\n", numel(mats));
fprintf("Options:\n");
fprintf("  1) Plot one file (choose by index)\n");
fprintf("  2) Make PDF with N animals per group\n");
choice = input("Enter 1 or 2 (default 2): ", "s");
if strlength(choice)==0, choice="2"; end

matFiles = fullfile({mats.folder}, {mats.name});

switch choice
    case "1"
        for i=1:numel(matFiles)
            fprintf("%3d. %s\n", i, mats(i).name);
        end
        idx = str2double(input("File number: ", "s"));
        if isnan(idx) || idx < 1 || idx > numel(matFiles)
            fprintf("Invalid index.\n");
            return;
        end
        plotOne(matFiles{idx}, [], false);
    otherwise
        nStr = input("How many per group? (default 3): ", "s");
        nPerGroup = str2double(nStr);
        if isnan(nPerGroup) || nPerGroup<=0, nPerGroup=3; end

        selected = pickNPerGroup(matFiles, nPerGroup, GROUP_PATTERNS);

        outPdf = fullfile(SIGNAL_DIR, sprintf("raw_traces_delta_sigma_%dpergroup.pdf", nPerGroup));
        fprintf("Saving: %s\n", outPdf);

        if exist(outPdf, "file"), delete(outPdf); end

        for i=1:numel(selected)
            fprintf("  plotting: %s\n", string(selected{i}));
            plotOne(selected{i}, outPdf, true);
        end

        fprintf("Done.\n");
end

% ===========================
% NESTED FUNCTIONS
% ===========================

    function plotOne(matPath, pdfPath, appendToPdf)
        [~, base, ~] = fileparts(matPath);   % base without extension
        base = string(base);
        csvPath = fullfile(SCORES_DIR, base + "_scored_scores_1Hz.csv");
        if ~isfile(csvPath)
            fprintf("⚠️  No scores CSV for %s, skipping.\n", base);
            return;
        end

        group = inferGroup(base, GROUP_PATTERNS);

        [eeg, emg, fs] = loadMatSignal(matPath);
        [time_s, states] = loadScores1Hz(csvPath);

        total_sec = numel(eeg) / fs;
        valid = time_s < total_sec;
        time_s = time_s(valid);
        states = states(valid);

        segs = compressStateSegments(time_s, states);
        warmIntervals = loadWarmingIntervals(WARMING_XLSX, base);

        [tD, delta_db] = bandpower1HzHilbert(eeg, fs, DELTA_BAND, POWER_BIN_SEC, POWER_SMOOTH_SEC);
        [tS, sigma_db] = bandpower1HzHilbert(eeg, fs, SIGMA_BAND, POWER_BIN_SEC, POWER_SMOOTH_SEC);

        % Downsample EEG for plotting
        n = numel(eeg);
        step = max(1, ceil(n / EEG_MAX_PLOT_POINTS));
        idx = 1:step:n;
        t_eeg = (idx-1) ./ fs;
        eeg_p = eeg(idx);

        % Figure + axes
        fig = figure('Color','w','Position',[100 100 1400 850]);
        tl = tiledlayout(fig, 4, 1, 'TileSpacing','compact', 'Padding','compact');

        ax1 = nexttile(tl,1); hold(ax1,'on');
        ax2 = nexttile(tl,2); hold(ax2,'on');
        ax3 = nexttile(tl,3); hold(ax3,'on');
        ax4 = nexttile(tl,4); hold(ax4,'on');

        % Plot traces
        plot(ax1, t_eeg, eeg_p, 'k', 'LineWidth', 0.4);
        ylabel(ax1, "EEG (\muV)");
        title(ax1, sprintf("%s  |  Group: %s  |  fs=%.1f Hz  |  duration=%.2f h", ...
            base, group, fs, total_sec/3600));

        if ~isempty(tD)
            plot(ax2, tD, delta_db, 'k', 'LineWidth', 1.0);
        end
        ylabel(ax2, "Delta power (dB)" + newline + "(0.5–5 Hz)");

        if ~isempty(tS)
            plot(ax3, tS, sigma_db, 'k', 'LineWidth', 1.0);
        end
        ylabel(ax3, "Sigma power (dB)" + newline + "(10–15 Hz)");

        % Warming strip panel
        ylim(ax4, [0 1]);
        yticks(ax4, []);
        xlabel(ax4, "Time (s)");
        ylabel(ax4, "Warming");
        % light grey background
        patch(ax4, [0 total_sec total_sec 0], [0 0 1 1], [0.95 0.95 0.95], ...
            'EdgeColor','none');

        % --- Shading: states on EEG/delta/sigma ---
        shadeSegments(ax1, segs, STATE_COLORS);
        shadeSegments(ax2, segs, STATE_COLORS);
        shadeSegments(ax3, segs, STATE_COLORS);

        % --- Shading: warming intervals across all panels ---
        shadeWarming(ax1, warmIntervals, [1.0 0.85 0.20], 0.20);
        shadeWarming(ax2, warmIntervals, [1.0 0.85 0.20], 0.20);
        shadeWarming(ax3, warmIntervals, [1.0 0.85 0.20], 0.20);
        shadeWarming(ax4, warmIntervals, [1.0 0.40 0.40], 0.60);

        % Optional "Warming ON" label on EEG axis
        for k=1:size(warmIntervals,1)
            ws = warmIntervals(k,1);
            text(ax1, ws, 0.98, "Warming ON", 'Units','normalized', ...
                'VerticalAlignment','top', 'HorizontalAlignment','left', 'FontSize',9);
        end

        % Cosmetics
        linkaxes([ax1 ax2 ax3 ax4], 'x');
        set([ax1 ax2 ax3 ax4], 'Box','off');

        % Legend (only show states present)
        uniqueStates = unique(states(:))';
        handles = gobjects(0); labels = strings(0);

        for st = sort(uniqueStates)
            if isKey(STATE_COLORS, double(st))
                rgba = STATE_COLORS(double(st));
                h = patch(ax1, nan, nan, rgba(1:3), 'EdgeColor','none', 'FaceAlpha', rgba(4));
                handles(end+1) = h; %#ok<AGROW>
                if isKey(STATE_NAMES, double(st))
                    labels(end+1) = string(STATE_NAMES(double(st)));
                else
                    labels(end+1) = string(st);
                end
            end
        end

        hWarm = patch(ax1, nan, nan, [1.0 0.85 0.20], 'EdgeColor','none', 'FaceAlpha', 0.20);
        handles(end+1) = hWarm;
        labels(end+1) = "Ambient warming window";

        legend(ax1, handles, labels, 'Location','northeast', 'FontSize',9);

        drawnow;

        % Save to PDF (append)
        if ~isempty(pdfPath)
            saveFigToPDF(fig, pdfPath, appendToPdf);
            close(fig);
        end
    end

    function name = matsName(p)
        [~, name, ext] = fileparts(p);
        name = name + ext;
    end
end

% ===========================
% HELPERS (separate functions)
% ===========================

function [eeg, emg, fs] = loadMatSignal(matPath)
S = load(matPath);
if ~isfield(S,'eeg') || ~isfield(S,'emg') || ~isfield(S,'eeg_frequency')
    error("MAT file %s must contain variables: eeg, emg, eeg_frequency", matPath);
end
eeg = double(S.eeg(:));
emg = double(S.emg(:));
fs  = double(S.eeg_frequency);
fs  = fs(1);
end

function [time_s, states] = loadScores1Hz(csvPath)
T = readtable(csvPath);

% numeric columns
isNum = varfun(@isnumeric, T, "OutputFormat","uniform");
numNames = T.Properties.VariableNames(isNum);
if isempty(numNames)
    error("No numeric columns found in %s", csvPath);
end

% pick score column: contains "score" or "state" else smallest nunique
scoreCol = "";
for i=1:numel(numNames)
    n = lower(string(numNames{i}));
    if contains(n,"score") || contains(n,"state")
        scoreCol = string(numNames{i});
        break;
    end
end
if scoreCol == ""
    % smallest unique count
    ucounts = zeros(1,numel(numNames));
    for i=1:numel(numNames)
        ucounts(i) = numel(unique(T.(numNames{i})));
    end
    [~,ix] = min(ucounts);
    scoreCol = string(numNames{ix});
end

states = int32(T.(scoreCol));
states = states(:);

% optional time column
timeCol = "";
for i=1:numel(numNames)
    if string(numNames{i}) == scoreCol, continue; end
    n = lower(string(numNames{i}));
    if contains(n,"time") || contains(n,"sec") || contains(n,"second")
        timeCol = string(numNames{i});
        break;
    end
end

if timeCol ~= ""
    time_s = double(T.(timeCol));
    time_s = time_s(:);
else
    time_s = (0:numel(states)-1)';
end

% ensure monotonic increasing; otherwise fallback to index time
if any(diff(time_s) <= 0)
    time_s = (0:numel(states)-1)';
end
end

function intervals = loadWarmingIntervals(xlsxPath, base)
intervals = zeros(0,2);
if strlength(xlsxPath)==0 || ~isfile(xlsxPath)
    return;
end

T = readtable(xlsxPath);

% find start/end columns
cols = string(T.Properties.VariableNames);
lc = lower(cols);
startIdx = find(contains(lc,"start"), 1, "first");
endIdx   = find(contains(lc,"finish") | (contains(lc,"end") & ~contains(lc,"start")), 1, "first");

if isempty(startIdx) || isempty(endIdx)
    fprintf("⚠️  XLSX missing start/end columns: %s\n", strjoin(cols, ", "));
    return;
end

% extract mouse number from filename: "mouse1" etc
tok = regexp(base, "mouse\s*(\d+)", "tokens", "once", "ignorecase");
if isempty(tok)
    fprintf("⚠️  Could not extract mouse number from %s\n", base);
    return;
end
mouseNum = tok{1};

% first column: mouse identifier; match number as a whole token
firstCol = cols(1);
firstVals = string(T.(firstCol));
mask = ~ismissing(firstVals) & ~cellfun(@isempty, regexp(firstVals, "\<" + mouseNum + "\>", "once"));
sub = T(mask, :);

if height(sub)==0
    fprintf("⚠️  No warming intervals found for mouse %s in %s\n", mouseNum, base);
    return;
end

sCol = cols(startIdx);
eCol = cols(endIdx);

pairs = zeros(0,2);
for r=1:height(sub)
    s = toSeconds(sub.(sCol)(r));
    e = toSeconds(sub.(eCol)(r));
    if isempty(s) || isempty(e), continue; end
    if e < s, tmp=s; s=e; e=tmp; end
    pairs(end+1,:) = [double(s) double(e)]; %#ok<AGROW>
end

if isempty(pairs), return; end

% de-dup + sort
pairs = unique(round(pairs,6), 'rows');
pairs = sortrows(pairs,1);
intervals = pairs;
end

function s = toSeconds(x)
% Accept numeric seconds, duration, datetime (time-of-day), or strings "hh:mm:ss" / "mm:ss" / numeric
s = [];
if isempty(x) || (isnumeric(x) && any(isnan(x)))
    return;
end

if isnumeric(x)
    s = double(x);
    return;
end

if isduration(x)
    s = seconds(x);
    return;
end

if isdatetime(x)
    % interpret as time-of-day
    s = hour(x)*3600 + minute(x)*60 + second(x);
    return;
end

% string-ish
str = strtrim(string(x));
if strlength(str)==0 || str=="NaN"
    return;
end

m = regexp(str, "^(?<hh>\d+):(?<mm>\d+):(?<ss>\d+(\.\d+)?)$", "names", "once");
if ~isempty(m)
    s = str2double(m.hh)*3600 + str2double(m.mm)*60 + str2double(m.ss);
    return;
end

m = regexp(str, "^(?<mm>\d+):(?<ss>\d+(\.\d+)?)$", "names", "once");
if ~isempty(m)
    s = str2double(m.mm)*60 + str2double(m.ss);
    return;
end

v = str2double(str);
if ~isnan(v)
    s = v;
end
end

function [t, mean_db] = bandpower1HzHilbert(eeg, fs, band, bin_sec, smooth_sec)
nyq = 0.5 * fs;
lo = max(0.001, band(1));
hi = min(band(2), nyq*0.999);

[b,a] = butter(4, [lo hi]/nyq, 'bandpass');
x = filtfilt(b, a, double(eeg(:)));

amp = abs(hilbert(x));
p = amp.^2;

total_sec = numel(p) / fs;
n_bins = floor(total_sec / bin_sec);
if n_bins <= 1
    t = []; mean_db = [];
    return;
end

edges_s = (0:n_bins) * bin_sec;
edges_i = floor(edges_s * fs) + 1;  % 1-based
edges_i(end) = min(edges_i(end), numel(p)+1);

% enforce strictly increasing + at least 2
edges_i = unique(edges_i, 'stable');
if numel(edges_i) < 2
    t = []; mean_db = [];
    return;
end

nb = numel(edges_i)-1;
mean_p = zeros(nb,1);

for k=1:nb
    i0 = edges_i(k);
    i1 = edges_i(k+1) - 1;
    if i0 > numel(p), break; end
    i1 = min(i1, numel(p));
    if i1 < i0
        mean_p(k) = NaN;
    else
        mean_p(k) = mean(p(i0:i1));
    end
end

t = (edges_i(1:end-1) - 1) ./ fs;         % seconds (0-based)
mean_db = 10 .* log10(mean_p + 1e-20);    % dB

if smooth_sec > 0
    sigma = smooth_sec / bin_sec; % in bins
    mean_db = gaussianSmooth1D(mean_db, sigma);
end
end

function y = gaussianSmooth1D(x, sigma)
if sigma <= 0
    y = x; return;
end
x = x(:);
m = max(1, round(3*sigma));
k = (-m:m)';
g = exp(-0.5*(k./sigma).^2);
g = g / sum(g);
y = conv(x, g, 'same');
end

function segs = compressStateSegments(time_s, states)
time_s = double(time_s(:));
states = double(states(:));
if isempty(time_s)
    segs = zeros(0,3);
    return;
end

segs = zeros(0,3);
s0 = states(1);
t0 = time_s(1);

for i=2:numel(states)
    if states(i) ~= s0
        t1 = time_s(i);
        segs(end+1,:) = [t0 t1 s0]; %#ok<AGROW>
        s0 = states(i);
        t0 = time_s(i);
    end
end
segs(end+1,:) = [t0 time_s(end)+1.0 s0];
end

function shadeSegments(ax, segs, STATE_COLORS)
% Shade segments using giant y-range so it stays valid even if ylim changes
for k=1:size(segs,1)
    a = segs(k,1); b = segs(k,2); st = segs(k,3);
    if isKey(STATE_COLORS, st)
        rgba = STATE_COLORS(st);
        patch(ax, [a b b a], [-1e9 -1e9 1e9 1e9], rgba(1:3), ...
            'FaceAlpha', rgba(4), 'EdgeColor','none');
    end
end
% Put shading behind lines
ch = ax.Children;
if ~isempty(ch)
    uistack(ch(1:end), 'bottom');
end
end

function shadeWarming(ax, intervals, rgb, alpha)
for k=1:size(intervals,1)
    a = intervals(k,1); b = intervals(k,2);
    patch(ax, [a b b a], [-1e9 -1e9 1e9 1e9], rgb, ...
        'FaceAlpha', alpha, 'EdgeColor','none');
end
end

function group = inferGroup(base, GROUP_PATTERNS)
b = char(base);
if ~isempty(regexp(b, GROUP_PATTERNS.WT,  'once', 'ignorecase'))
    group = "WT";
elseif ~isempty(regexp(b, GROUP_PATTERNS.APP, 'once', 'ignorecase'))
    group = "APP";
else
    group = "UNK";
end
end

function selected = pickNPerGroup(matFiles, nPerGroup, GROUP_PATTERNS)
byGroup = struct('WT',{{}}, 'APP',{{}}, 'UNK',{{}});
for i=1:numel(matFiles)
    [~, nm] = fileparts(matFiles{i});
    g = inferGroup(string(nm), GROUP_PATTERNS);
    byGroup.(g){end+1} = matFiles{i}; %#ok<AGROW>
end

selected = {};
selected = [selected, byGroup.WT(1:min(nPerGroup, numel(byGroup.WT)))];
selected = [selected, byGroup.APP(1:min(nPerGroup, numel(byGroup.APP)))];
end

function saveFigToPDF(fig, pdfPath, append)
% R2020a+ supports exportgraphics(...,'Append',true)
try
    exportgraphics(fig, pdfPath, 'ContentType','vector', 'Append', append);
catch
    % fallback: print to temp then append via ghostscript-like not available
    % simplest fallback: overwrite (no append)
    warning("exportgraphics append not supported here. Saving single page only.");
    print(fig, pdfPath, '-dpdf', '-painters');
end
end
