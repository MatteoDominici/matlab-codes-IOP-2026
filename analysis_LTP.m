% =========================================================================
% This script processes a set of MEA recordings stored as .h5 (HDF5) files
% acquired during an LTP protocol. It automatically:
%
% 1) Loads two baseline control recordings (Control1 and Control2) and
%    performs spike detection and stimulation-event extraction using
%    AnalysisFile (band-pass filtering + threshold-based spike detection)
%
% 2) Computes, for each electrode/channel, the evoked response strength as
%    the area under the peristimulus time histogram (PSTH) within a fixed
%    post-stimulus window (winPSTH) and using a chosen binning (binSize).
%    The baseline evoked response (Area_BC) is defined as the mean of the
%    two control PSTH areas
%
% 3) Identifies post-LTP recordings by parsing filenames, keeps only the
%    expected time points (e.g., 0–180 min), and repeats the same analysis
%    to obtain PSTH areas at each post-LTP time point
%
% 4) Computes the per-channel efficacy change (DeltaMatrix) at each time
%    point as a relative variation from baseline. It and also computes a 
%    baseline stability index S from Control1 vs Control2 to quantify 
%    recording reliability per channel
%
% 5) Classifies electrodes into “red” (potentiated) or “blue” (non-potentiated)
%    based on baseline activity (Area_BC), baseline stability (S), and mean
%    post-LTP delta exceeding a stability-scaled threshold
%
% 6) Computes an additional pre-point (Δ_pre, corresponding to -37 min) from
%    Control2 vs baseline, prepends it to the delta matrices, and updates the
%    time labels so the saved results are self-consistent
%
% 7) Saves all key variables (baseline data, deltas, classifications, timing)
%    into a .mat file for later group-level aggregation across experiments
%
% 8) Computes a Potentiation Index (PI) over time, defined as the fraction of
%    red channels meeting the potentiation criterion at each time point over
%    the total number of classified channels (red + blue)
% 
% 9) Generates a MEA-layout plot showing Δ time courses for each
%    electrode (with stimulated electrodes highlighted from datasetName)
% =========================================================================

%% ========================================================================
%  ===== PARAMETERS =====

% Directory containing the recording files (.h5)
dataDir = 'C:\Users\matte\Documents\Università\Magistrale\Secondo anno\Secondo semestre\TESI LM\Hyperione\LTP\45188_17DIV';

fs = 20000; % sampling frequency (Hz)
binSize = 0.005; % PSTH bin width (s)
winPSTH = [0 0.250]; % PSTH analysis window relative to stimulation onset (s)
k = 7; % spike detection threshold multiplier
bpf = [300 3000]; % band-pass filter cutoff frequencies (Hz)

% Expected post-LTP time points (minutes)
post_times_expected = [0 30 60 90 120 180];

% Dataset identifier (used for labeling and for extracting stimulated channels)
datasetName = '45188_17DIV_(63-64)';

% Output directory for saving results and figures
outDir = 'C:\Users\matte\Documents\Università\Pubblicazione Articolo\Codici';
if ~exist(outDir,'dir'), mkdir(outDir); end

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
%  ===== ANALYSIS (Spike Detection, PSTH Area, Efficacy, Stability) =====

% Retrieve all HDF5 files in the data directory
files = dir(fullfile(dataDir, '*.h5'));
if isempty(files), error('No files .h5 in %s', dataDir); end

% Identify indices of the 2 baseline control recordings (Control 1 and 2)
c1_idx = find(contains({files.name}, 'Control1'));
c2_idx = find(contains({files.name}, 'Control2'));
if isempty(c1_idx) || isempty(c2_idx)
    error('No files "Control1/Control2" in %s', dataDir);
end

% Build full paths to the first occurrences of Control1 and Control2
fileControl1 = fullfile(files(c1_idx(1)).folder, files(c1_idx(1)).name);
fileControl2 = fullfile(files(c2_idx(1)).folder, files(c2_idx(1)).name);

% ------------------------------------------------------------------------
% Identify post-LTP files from the filename pattern and select expected times

post_files = {}; post_times = [];
for f = 1:numel(files)
    % Extract the post-LTP time (in minutes) from the filename
    tok = regexp(files(f).name, 'WF_0\.2Hz_(\d+)min', 'tokens','once');
    if ~isempty(tok)
        tmin = str2double(tok{1});
        % Retain only time points included in the expected post-LTP list
        if ismember(tmin, post_times_expected)
            post_files{end+1} = fullfile(files(f).folder, files(f).name);
            post_times(end+1) = tmin;
        end
    end
end

% Sort post-LTP files by time in ascending order
[post_times, idx] = sort(post_times);
post_files = post_files(idx);

% ------------------------------------------------------------------------
% Perform spike detection once for the two baseline control recordings
[spikes_c1, Ch_c1, tpOn_c1] = AnalysisFile(fileControl1, fs, k, bpf);
[spikes_c2, Ch_c2, tpOn_c2] = AnalysisFile(fileControl2, fs, k, bpf);

% ------------------------------------------------------------------------
% Compute PSTH areas for Control1 and Control2 using the chosen binning/window
Area_c1 = computeAreaPSTH(spikes_c1, tpOn_c1, binSize, winPSTH);
Area_c2 = computeAreaPSTH(spikes_c2, tpOn_c2, binSize, winPSTH);
Area_BC = (Area_c1 + Area_c2)/2; % baseline control area as mean

% ------------------------------------------------------------------------
% Perform spike detection for each post-LTP recording and store results
PostSpikes = struct; % temporary structure
for f = 1:numel(post_files)
    [spikes_post, Ch_post, tpOn_post] = AnalysisFile(post_files{f}, fs, k, bpf);
    PostSpikes(f).spikes   = spikes_post;
    PostSpikes(f).Channels = Ch_post;
    PostSpikes(f).tpOn     = tpOn_post;
end

% ------------------------------------------------------------------------
% Compute post-LTP PSTH areas and assemble a results structure per time point
Results = struct;
for f = 1:numel(PostSpikes)
    % Retrieve per-file content from the temporary structure
    spikes_post = PostSpikes(f).spikes;
    tpOn_post   = PostSpikes(f).tpOn;
    Ch_post     = PostSpikes(f).Channels;
    
    % Compute PSTH area for the current post-LTP recording
    Area_post = computeAreaPSTH(spikes_post, tpOn_post, binSize, winPSTH);
    
    % Store output fields for downstream processing and plotting
    Results(f).Area     = Area_post;
    Results(f).spikes   = spikes_post;
    Results(f).Channels = Ch_post;
    Results(f).tpOn     = tpOn_post;
    Results(f).time     = post_times(f);
end

% ------------------------------------------------------------------------
% DELTA (efficacy) computation

% Collect all post-LTP time stamps and create a sorted time index
allTimes = [Results.time];
[timeLabelsSorted, idxSort] = sort(allTimes);

% Number of post-LTP time points
nTimes   = numel(Results);

% Channel list derived from Control1 (reference channel order)
allCh    = Ch_c1(:)';

% Preallocate delta matrix: rows = channels, columns = time points
DeltaMatrix = nan(numel(allCh), nTimes);

% Compute delta for each channel at each time point
for t = 1:nTimes
    % Channel labels available in the current post-LTP file
    ch_t = Results(t).Channels(:)';
    for ic = 1:numel(allCh)
        % Current reference channel identifier
        ch = allCh(ic);
        % Find the channel position in the current post-LTP channel list
        pos = find(ch_t == ch, 1);
        % Assign NaN if channel is missing or baseline is zero (avoid division)
        if isempty(pos) || Area_BC(ic)==0
            DeltaMatrix(ic, t) = NaN;
        else
            % Relative change: (post - baseline) / baseline
            DeltaMatrix(ic, t) = (Results(t).Area(pos) - Area_BC(ic)) ./ Area_BC(ic);
        end
    end
end

% Baseline stability index between Control1 and Control2 (per channel)
S = abs((Area_c2(:) - Area_c1(:)) ./ (Area_c2(:) + Area_c1(:) + eps));

%% ========================================================================
%  ===== Red/Blue Electrodes Selection =====

% Initialize logical masks for electrode classification
isRed  = false(numel(allCh),1);
isBlue = false(numel(allCh),1);

for ic = 1:numel(allCh)
    % Post-LTP delta values for the current channel, ordered by time
    d = DeltaMatrix(ic, idxSort);
    
    % Baseline activity criterion (minimum baseline evoked activity)
    cond1 = Area_BC(ic) >= 0.5;

    % Baseline stability criterion between Control1 and Control2
    cond2 = S(ic) <= 0.2;

    % Compute mean post-LTP delta while ignoring missing values
    if all(isnan(d))
        meanDelta = NaN;
    else
        meanDelta = mean(d, 'omitnan');
    end

    % Potentiation criterion: mean post-LTP delta must exceed a stability-scaled threshold
    cond3 = ~isnan(meanDelta) && (meanDelta >= 1.2 * S(ic));

    % Assign channel class based on the combined criteria
    if cond1 && cond2 && cond3
        isRed(ic) = true;
    elseif cond1 && cond2 && ~cond3
        isBlue(ic) = true;
    end
end

% Extract delta matrices restricted to red and blue channel subsets
Delta_red_all  = DeltaMatrix(isRed, :);
Delta_blue_all = DeltaMatrix(isBlue, :);

%% ========================================================================
%  ===== Compute Δ_pre and include it in the saved outputs =====

% Compute the "pre" delta for each channel.
% Here Δ_pre quantifies the relative change between Control2 and the baseline
% control average (Area_BC): (Control2 - Baseline) / Baseline.
Delta_pre = (Area_c2(:) - Area_BC(:)) ./ Area_BC(:);

% Prepend the Δ_pre column to the red/blue delta matrices so that:
%   - column 1 = pre (-37 min)
%   - columns 2..end = post-LTP time points (0..180 min)
Delta_red_all  = [Delta_pre(isRed),  DeltaMatrix(isRed, :)];
Delta_blue_all = [Delta_pre(isBlue), DeltaMatrix(isBlue, :)];

% Update time labels accordingly so the saved .mat is self-consistent
timeLabelsSorted = [-37, timeLabelsSorted];

%% ========================================================================
%  ===== Saving =====

save(fullfile(outDir,'45188_17DIV.mat'), ...
     'Area_c1','spikes_c1','Ch_c1','tpOn_c1', ...
     'Area_c2','spikes_c2','Ch_c2','tpOn_c2', ...
     'Area_BC','Results','DeltaMatrix','layout','binSize','winPSTH', ...
     'post_times','S','timeLabelsSorted', ...
     'Delta_red_all','Delta_blue_all');

%% ========================================================================
%  ===== Potentiation Index (PI) at different time points =====

% Compute the Potentiation Index (PI) at each time point.
% PI is defined as the fraction of "valid potentiated" red channels
% (Δ >= 1.2*S) over the total number of classified channels (red + blue).

nTimes = numel(timeLabelsSorted);
PI = nan(1, nTimes);

for t = 1:nTimes
    % Delta values for red channels at time t
    Delta_red_t = DeltaMatrix(isRed, t);
    S_red = S(isRed);
    
    % Potentiation criterion: Δ >= 1.2 * S
    valid_red = ~isnan(Delta_red_t) & (Delta_red_t >= 1.2 .* S_red);
    n_valid = sum(valid_red);

    % Total number of classified channels (red + blue)
    n_total = sum(isRed) + sum(isBlue);

    % Compute PI for the current time point
    PI(t) = n_valid / n_total;
end

% Print results to the command window
fprintf('\n=== POTENTIATION INDEX (PI) ===\n');
for t = 1:nTimes
    fprintf('Time %3d min -> PI = %.3f\n', timeLabelsSorted(t), PI(t));
end
fprintf('===============================\n\n');

%% ========================================================================
%  ===== Aggregated ΔA over time (plot in MEA layout) =====
figure('Name','Delta vs Time (MEA layout)','Units','normalized','OuterPosition',[0 0 1 1]);
tiledlayout(8,8,'TileSpacing','compact','Padding','compact');

% Extract stimulated electrode IDs from datasetName (e.g., "(63-64)")
stimCh = regexp(datasetName, '\((\d+)-(\d+)\)', 'tokens');
if ~isempty(stimCh)
    stimCh = str2double(stimCh{1});
else
    stimCh = [];
end

% Add the "pre" point to the time axis and to the delta vectors
% Δ_pre is computed consistently with the saved version above
Delta_pre = (Area_c2(:) - Area_BC(:)) ./ Area_BC(:);
post_times_ext = [-37, post_times];

for row = 1:8
    for col = 1:8
        chNumber = layout(row,col);
        nexttile;

        % Skip empty positions in the MEA layout
        if isnan(chNumber), axis off; continue; end

        % Find channel index in the reference channel list
        ic = find(allCh==chNumber,1);

        % If the channel is not present, show only the label (optional highlight)
        if isempty(ic)
            axis off;
            if ismember(chNumber, stimCh)
                title(sprintf('Ch %d',chNumber),'FontSize',8,'BackgroundColor',[1 1 0]);
            else
                title(sprintf('Ch %d',chNumber),'FontSize',8);
            end
            continue;
        end

        % Channel delta time course (sorted by time) + prepend pre value
        d = DeltaMatrix(ic, idxSort);
        d_all = [Delta_pre(ic), d];

        % Choose line color based on electrode classification
        if isRed(ic)
            colLine = [1 0 0];
        elseif isBlue(ic)
            colLine = [0 0 1];
        else
            colLine = [0.6 0.6 0.6];
        end
        
        % Plot Δ over time for the current electrode
        plot(post_times_ext, d_all, '-o', 'LineWidth', 1.2, 'Color', colLine); hold on;

        % Reference lines: 0 (no change) and S (stability index for that channel)
        yline(0, ':k');
        yline(S(ic), ':', 'Color', [0 0.4 0], 'LineWidth', 1.2);

        % Title (highlight the 2 stimulated electrodes)
        if ismember(chNumber, stimCh)
            title(sprintf('Electrode %d', chNumber),'FontSize',8,'BackgroundColor',[1 1 0]);
        else
            title(sprintf('Electrode %d', chNumber),'FontSize',8);
        end

        % Axis formatting
        xlim([min(post_times_ext) max(post_times_ext)]);
        xticks(post_times_ext);
        xticklabels(arrayfun(@num2str, post_times_ext, 'UniformOutput', false));
        set(gca,'XTickLabelRotation',0);
        set(gca,'FontSize',7,'Box','off'); 
        grid on;

        % Axis labels
        xlabel('Time [min]','FontSize',7);
        ylabel('\DeltaA','FontSize',7,'Interpreter','tex');

        % Automatic y-limits based on the plotted data
        if all(isnan(d_all))
            ylim([-1 1]);
        else
            Ymin = min(d_all,[],'omitnan'); 
            Ymax = max(d_all,[],'omitnan');
            if Ymax==Ymin
                ylim([Ymin-0.5, Ymax+0.5]);
            else
                ylim([Ymin-0.1*(Ymax-Ymin), Ymax+0.1*(Ymax-Ymin)]);
            end
        end
    end
end

% Global title and export
safeName = strrep(datasetName, '_', '\_');
sgtitle(['$\Delta A$ evoked response - ' safeName], 'Interpreter', 'latex');
exportgraphics(gcf, fullfile(outDir,'Delta_vsTime_45188_17DIV.jpg'), 'Resolution',300);
close(gcf);

% =========================================================================
%% ========================================================================
%  ===== Local function for PSTH area =====

function Area = computeAreaPSTH(spike_times_sec, tpOn, binSize, winPSTH)
    % Compute the PSTH area (integral of the mean firing rate) in a given window
    %
    % For each channel:
    % 1) align spikes to each stimulation onset (tpOn)
    % 2) bin spikes within the analysis window (winPSTH) using binSize
    % 3) normalize to obtain firing rate in Hz
    % 4) integrate over time to get the PSTH area

    nCh = length(spike_times_sec);
    Area = zeros(1,nCh);

    for ch = 1:nCh
        sp = spike_times_sec{ch};

        % If no spikes are present, area is set to 0
        if isempty(sp)
            Area(ch) = 0;
            continue;
        end
        
        % Define PSTH bin edges within the selected window
        edges = winPSTH(1):binSize:winPSTH(2);
        nBins = length(edges)-1;

        % Accumulate spike counts across stimuli
        psth_counts = zeros(1, nBins);

        for i = 1:length(tpOn)
            stim_time = tpOn(i);

            % Collect spikes occurring in the window relative to the stimulus
            spikes_in_window = sp(sp >= stim_time + winPSTH(1) & ...
                                  sp <  stim_time + winPSTH(2)) - stim_time;
            if isempty(spikes_in_window), continue; end

            % Bin spikes for the current stimulus and add to the accumulator
            counts = histcounts(spikes_in_window, edges);
            psth_counts = psth_counts + counts;
        end

        % Normalize: average across stimuli and convert counts/bin to Hz
        psth_norm = psth_counts / length(tpOn) / binSize;

        % Area under the PSTH curve (integral over time)
        Area(ch) = sum(psth_norm) * binSize;
    end
end
