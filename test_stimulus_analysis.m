%% ========================================================================
%  MEA SPIKE PIPELINE + TR/TB RATE ANALYSIS + SIGMOID FIT
%
%  This script performs a 3-stage analysis on MEA recordings:
%   1) Spike extraction + preprocessing per file (via AnalysisFile)
%   2) Spike-rate computation in two windows (Tr and Tb) + 8x8 diagnostic maps
%   3) Sigmoid fitting (k, t0, Delta) on selected MEA channels for Tr and Tb
%
%  The workflow is organized as 3 loops over the input files. Intermediate
%  results are stored in a struct ('results') indexed by filename.
% =========================================================================

%% ========================================================================
%  ===== PARAMETERS =====

% Directory containing the recording files (.h5) of the specified dataset
dataDir = 'E:\Matteo_tesiLM\frequency_WF_baseline_experiments\40628_11DIV';

fs = 20000; % sampling frequency (Hz)
binSize = 0.005; % PSTH bin width (s)
winPSTH = [0 0.250]; % PSTH analysis window relative to stimulation onset (s)
k = 7; % spike detection threshold multiplier
bpf = [300 3000]; % band-pass filter cutoff frequencies (Hz)

% Discover data files (.h5) matching the used naming convention (WF_fr)
files = dir(fullfile(dataDir, '*WF_fr(*.h5'));
if isempty(files)
    error('No files found.');
end

%% ========================================================================
%  ===== MEA LAYOUT (8x8) =====

% Electrode layout matrix (NaN indicates no physical electrode position)
layout = [
     NaN, 21, 31, 41, 51, 61, 71, NaN;
     12, 22, 32, 42, 52, 62, 72, 82;
     13, 23, 33, 43, 53, 63, 73, 83;
     14, 24, 34, 44, 54, 64, 74, 84;
     15, 25, 35, 45, 55, 65, 75, 85;
     16, 26, 36, 46, 56, 66, 76, 86;
     17, 27, 37, 47, 57, 67, 77, 87;
     NaN, 28, 38, 48, 58, 68, 78, NaN
];

% Exclude reference electrode from downstream analysis
layout(layout == 15) = NaN;

%% ========================================================================
%  ===== 1st LOOP: Spike detection (via AnalysisFile) =====
%
%  For each recording file, run the extraction + filtering + spike detection
%  pipeline. Store:
%   - spike_times_sec : cell array (one cell per channel)
%   - Channels        : channel list/IDs
%   - tpOn            : stimulation onset timestamps (s); empty if no stimuli

results = struct;

for f = 1:numel(files)

    filename = fullfile(files(f).folder, files(f).name);
    fprintf('\n=== File %s ===\n', files(f).name);

    % Run the full preprocessing and spike detection pipeline.
    % OUTPUTS are assumed to be:
    %   spike_times_sec : {nChannels x 1} cell, spike timestamps in seconds
    %   Channels        : vector of MEA channel IDs
    %   tpOn            : vector of stimulus onset times (s)
    [spike_times_sec, Channels, tpOn] = AnalysisFile(filename, fs, k, bpf);

    % Store results using a valid MATLAB field name derived from filename
    varName = matlab.lang.makeValidName(files(f).name);
    results.(varName).filename        = files(f).name;
    results.(varName).Channels        = Channels;
    results.(varName).spike_times_sec = spike_times_sec;
    results.(varName).tpOn            = tpOn;  % empty if no stimuli are present
end

%% ========================================================================
%  ===== 2nd LOOP: Tr / Tb spike-rate analysis + 8x8 plots =====
%
%  For each file:
%   - Estimate stimulation period T from tpOn
%   - Compute per-pulse spike-rate (Hz) in:
%       Tr : [tpOn(n), tpOn(n)+Tr]
%       Tb : [tpOn(n)+Tr, tpOn(n+1)]
%   - Smooth per-channel SR using an adaptive moving mean window
%   - Estimate "decay time" as the earliest time after which the smoothed
%     spike-rate stays below a fixed threshold for the remaining pulses
%   - Produce and save 8x8 tiled maps (raw + smoothed + decay marker)

outDir = 'E:\Matteo_tesiLM\frequency_WF_baseline_experiments\40628_11DIV';
if ~exist(outDir,'dir'), mkdir(outDir); end

decayThreshold = 0.5; % threshold used to define decay time (Hz)

fileNames = fieldnames(results);

for f = 1:numel(fileNames)

    res = results.(fileNames{f});

    spike_times_sec = res.spike_times_sec;
    Channels        = res.Channels;
    tpOn            = res.tpOn;

    % If no stimuli are present, skip SR analysis
    if isempty(tpOn) || numel(tpOn) < 2
        warning('No valid tpOn found for %s. Skipping SR analysis.', res.filename);
        continue;
    end

    % Estimate stimulation period T from onset differences
    T_array = diff(tpOn);
    T = median(T_array); % mean stimulation period (s)
    nStimuli = length(tpOn) - 1; % number of inter-stimulus intervals

    % Save period in results for later use
    results.(fileNames{f}).T = T;

    % ---------------------------------------------------------------------
    % Build mapping: MEA channel ID -> index in the Channels array
    % This is used to place data correctly in the 8x8 layout plots
    channelToIndex = containers.Map('KeyType','double','ValueType','double');
    for ch = 1:length(Channels)
        if ~isnan(Channels(ch))
            channelToIndex(Channels(ch)) = ch;
        end
    end

    %% ============================= Tr analysis ==========================
    % Time base for Tr is aligned to stimulus onset (tpOn(n)).
    timeVector_Tr = tpOn(1:nStimuli);

    SR_Tr = zeros(length(Channels), nStimuli); % raw SR per pulse
    SRs_Tr = zeros(length(Channels), nStimuli); % smoothed SR per pulse
    decayTimes_Tr = NaN(1, length(Channels)); % decay time estimate

    for ch = 1:length(Channels)

        spikes = spike_times_sec{ch};
        if isempty(spikes), continue; end

        % Compute raw SR in window [tpOn(n), tpOn(n)+Tr]
        for n = 1:nStimuli
            t_start = tpOn(n);
            t_end   = t_start + Tr;

            s_in = spikes(spikes >= t_start & spikes < t_end);
            SR_Tr(ch, n) = numel(s_in) / Tr;
        end

        raw_sr = SR_Tr(ch,:);
        if all(raw_sr==0), continue; end

        % Adaptive smoothing window:
        % - base_window uses average derivative magnitude
        % - decay_based uses total duration / (5*T)
        % Final window is constrained to [3, 20]
        d_sr        = abs(diff(raw_sr));
        avg_deriv   = mean(d_sr,'omitnan');
        base_window = round(avg_deriv / 0.05);
        decay_based = round((timeVector_Tr(end)-timeVector_Tr(1))/(5*T));
        smoothWindowSize = min(20, max(3, round((base_window+decay_based)/2)));

        smooth_sr    = movmean(raw_sr, smoothWindowSize);
        SRs_Tr(ch,:) = smooth_sr;

        % Decay time definition:
        % first timepoint n such that all remaining smoothed SR values
        % stay below decayThreshold.
        for n = 1:length(smooth_sr)
            if all(smooth_sr(n:end) < decayThreshold)
                decayTimes_Tr(ch) = timeVector_Tr(n);
                break;
            end
        end
    end

    % ---------------------- Plot Tr (8x8 tiled map) ----------------------
    figure('Units','normalized','OuterPosition',[0 0 1 1]);
    tiledlayout(8,8,'TileSpacing','compact','Padding','compact');

    pulseIdx_Tr = 1:nStimuli;

    for row = 1:8
        for col = 1:8

            chNumber = layout(row,col);
            nexttile;

            % Empty MEA positions: no axes
            if isnan(chNumber)
                axis off;
                continue;
            end

            % Channel not found in current file: mark as blank tile
            if ~isKey(channelToIndex, chNumber)
                fill([0 1 1 0],[0 0 1 1],[1 1 1]); axis off;
                title(sprintf('Ch %d', chNumber), 'FontSize', 7);
                continue;
            end

            idx         = channelToIndex(chNumber);
            raw_sr      = SR_Tr(idx,:);
            smoothed_sr = SRs_Tr(idx,:);

            % Dynamic y-limits to improve readability across channels
            yMin   = min([raw_sr, smoothed_sr]);
            yMax   = max([raw_sr, smoothed_sr]);
            yRange = yMax - yMin;
            if yRange == 0, yRange = 1; end
            yLims  = [max(0, yMin - 0.15*yRange), yMax + 0.15*yRange];

            % Raw + smoothed
            plot(pulseIdx_Tr, raw_sr, 'x', 'Color',[0.4 0.4 0.4], 'MarkerSize', 4); hold on;
            plot(pulseIdx_Tr, smoothed_sr, 'b-', 'LineWidth', 1.2);

            % Decay marker (vertical line at decay pulse index)
            if ~isnan(decayTimes_Tr(idx))
                [~, decayIdx] = min(abs(timeVector_Tr - decayTimes_Tr(idx)));
                xline(decayIdx, 'r--', 'LineWidth', 1);
            end

            xlabel('Pulse #', 'FontSize', 6);
            ylabel('Rate (Hz)', 'FontSize', 6);
            title(sprintf('Ch %d', chNumber), 'FontSize', 7);
            set(gca,'FontSize',6);
            ylim(yLims);
        end
    end

    % Clean base file name for titles and output filenames
    [~, fname, ~] = fileparts(res.filename);
    fname = regexprep(fname,'__D_\d+$','');

    sgtitle(sprintf('Spike rate (raw + smoothed) + decay — Tr — %s', fname), ...
        'Interpreter','none');

    print(gcf, fullfile(outDir, [fname '_Tr.jpg']), '-djpeg', '-r300');
    close;

    %% ============================= Tb analysis ==========================
    % Tb window starts immediately after Tr and ends at next stimulus onset
    timeVector_Tb = tpOn(1:nStimuli) + Tr;

    SR_Tb        = zeros(length(Channels), nStimuli);
    SRs_Tb       = zeros(length(Channels), nStimuli);
    decayTimes_Tb = NaN(1, length(Channels));

    for ch = 1:length(Channels)

        spikes = spike_times_sec{ch};
        if isempty(spikes), continue; end

        % Compute raw SR in window [tpOn(n)+Tr, tpOn(n+1)]
        for n = 1:nStimuli
            t_start = tpOn(n) + Tr;
            t_end   = tpOn(n+1);

            % Window duration may vary slightly: ensure strictly positive
            winDur = max(eps, t_end - t_start);

            s_in = spikes(spikes >= t_start & spikes < t_end);
            SR_Tb(ch, n) = numel(s_in) / winDur;
        end

        raw_sr = SR_Tb(ch,:);
        if all(raw_sr==0), continue; end

        % Adaptive smoothing window
        d_sr        = abs(diff(raw_sr));
        avg_deriv   = mean(d_sr,'omitnan');
        base_window = round(avg_deriv / 0.05);
        decay_based = round((timeVector_Tb(end)-timeVector_Tb(1))/(5*T));
        smoothWindowSize = min(20, max(3, round((base_window+decay_based)/2)));

        smooth_sr    = movmean(raw_sr, smoothWindowSize);
        SRs_Tb(ch,:) = smooth_sr;

        % Decay time definition (same criterion as Tr)
        for n = 1:length(smooth_sr)
            if all(smooth_sr(n:end) < decayThreshold)
                decayTimes_Tb(ch) = timeVector_Tb(n);
                break;
            end
        end
    end

    % ---------------------- Plot Tb (8x8 tiled map) ----------------------
    figure('Units','normalized','OuterPosition',[0 0 1 1]);
    tiledlayout(8,8,'TileSpacing','compact','Padding','compact');

    pulseIdx_Tb = 1:nStimuli;

    for row = 1:8
        for col = 1:8

            chNumber = layout(row,col);
            nexttile;

            if isnan(chNumber)
                axis off;
                continue;
            end

            if ~isKey(channelToIndex, chNumber)
                fill([0 1 1 0],[0 0 1 1],[1 1 1]); axis off;
                title(sprintf('Ch %d', chNumber),'FontSize',7);
                continue;
            end

            idx         = channelToIndex(chNumber);
            raw_sr      = SR_Tb(idx,:);
            smoothed_sr = SRs_Tb(idx,:);

            yMin   = min([raw_sr, smoothed_sr]);
            yMax   = max([raw_sr, smoothed_sr]);
            yRange = yMax - yMin;
            if yRange == 0, yRange = 1; end
            yLims  = [max(0, yMin - 0.15*yRange), yMax + 0.15*yRange];

            plot(pulseIdx_Tb, raw_sr, 'x', 'Color',[0.4 0.4 0.4], 'MarkerSize', 4); hold on;
            plot(pulseIdx_Tb, smoothed_sr, 'b-', 'LineWidth', 1.2);

            if ~isnan(decayTimes_Tb(idx))
                [~, decayIdx] = min(abs(timeVector_Tb - decayTimes_Tb(idx)));
                xline(decayIdx, 'r--', 'LineWidth', 1);
            end

            xlabel('Pulse #', 'FontSize', 6);
            ylabel('Rate (Hz)', 'FontSize', 6);
            title(sprintf('Ch %d', chNumber),'FontSize',7);
            set(gca,'FontSize',6);
            ylim(yLims);
        end
    end

    sgtitle(sprintf('Spike rate (raw + smoothed) + decay — Tb — %s', fname), ...
        'Interpreter','none');

    print(gcf, fullfile(outDir, [fname '_Tb.jpg']), '-djpeg', '-r300');
    close;

    % ---------------------------------------------------------------------
    % Save SR matrices for subsequent sigmoid fitting (3rd loop)
    results.(fileNames{f}).SR_Tr  = SR_Tr;
    results.(fileNames{f}).SRs_Tr = SRs_Tr;
    results.(fileNames{f}).SR_Tb  = SR_Tb;
    results.(fileNames{f}).SRs_Tb = SRs_Tb;
end

%% ========================================================================
%  ===== Selected channels for fitting (MEA channel IDs) =====

% Only these channels will be passed to the sigmoid fitting routine
% These are those electrodes exhibithing a spike rate drop at 0.5 or 1.0 Hz
% and above
selectedChannels_Tr = [14 17 23 24 25 26 27 35 42 53];
selectedChannels_Tb = [17 24 25 26 35 53];

%% ========================================================================
%  ===== 3rd LOOP: Sigmoid fit (k, t0, Delta) for Tr and Tb =====
%
%  For each file (sorted by stimulation frequency):
%   - Extract frequency from filename (fr(XXHz))
%   - Compute T = 1/f
%   - Fit sigmoid to SR time courses for the selected channels, separately
%     for Tr and Tb, using fitSigmoidForSelectedChannels()
%   - Save fitResults struct to disk

% -------------------------------------------------------------------------
% Sort files by stimulation frequency extracted from the results field name
fileNames = fieldnames(results);

freqVals = nan(1, numel(fileNames));
for i = 1:numel(fileNames)
    tok = regexp(fileNames{i}, 'fr\((\d+\.?\d*)Hz\)', 'tokens', 'once');
    if ~isempty(tok)
        freqVals(i) = str2double(tok{1});
    end
end
[~, sortIdx] = sort(freqVals, 'descend');
fileNames = fileNames(sortIdx);

% Silence common fitting warnings
warning('off','MATLAB:singularMatrix');
warning('off','optimlib:lsqncommon:FiniteDiffStepSize');
warning('off','stats:nlparci:JacobianNotFullRank');

% -------------------------------------------------------------------------
for f = 1:numel(fileNames)

    res = results.(fileNames{f});

    % Clean base name (used for plots and struct field names)
    fname     = regexprep(res.filename, '__D_\d+$', '');
    validName = matlab.lang.makeValidName(fname);

    % Extract stimulation frequency (Hz) from filename
    tok = regexp(fname, 'fr\((\d+\.?\d*)Hz\)', 'tokens', 'once');
    if isempty(tok)
        fprintf('No frequency found in %s\n', fname);
        continue;
    end
    freqVal = str2double(tok{1});
    T       = 1 / freqVal;     % stimulation period (s)
    fitResults.(validName).T = T;

    fprintf('\n=== File %s (freq=%.2f Hz, T=%.3f s) ===\n', ...
        fname, freqVal, T);

    % Ensure SR fields exist (generated in the 2nd loop)
    if ~isfield(res,'SRs_Tr') || ~isfield(res,'SR_Tr') || ...
       ~isfield(res,'SRs_Tb') || ~isfield(res,'SR_Tb')
        warning('Missing SR fields for %s. Run the 2nd loop first.', fname);
        continue;
    end

    %% ---------------------------- FIT Tr ----------------------------
    fprintf('--- FIT Tr ---\n');
    fitResults.(validName).Tr = fitSigmoidForSelectedChannels( ...
        res.Channels, ...
        res.SR_Tr, ...
        res.SRs_Tr, ...
        selectedChannels_Tr, ...
        T, outDir, fname, 'Tr');

    %% ---------------------------- FIT Tb ----------------------------
    fprintf('--- FIT Tb ---\n');
    fitResults.(validName).Tb = fitSigmoidForSelectedChannels( ...
        res.Channels, ...
        res.SR_Tb, ...
        res.SRs_Tb, ...
        selectedChannels_Tb, ...
        T, outDir, fname, 'Tb');

end

% -------------------------------------------------------------------------
% Save final fit results to disk
save(fullfile(outDir,'fitResults_40628_11DIV_k_t0_Delta.mat'), 'fitResults');
