% =========================================================================
% This script performs a post-processing and aggregation analysis of sigmoid
% fit results obtained from single-electrode LTP experiments. Specifically:
% 1) Loads multiple .mat files containing previously computed sigmoid-fit
%    parameters (fitResults structures) for Tr and Tb response windows.
% 2) Groups datasets by stimulation frequency and separates Tr and Tb data.
% 3) Applies quality-control criteria to exclude unreliable fits
% 4) Removes duplicate entries and performs limited manual curation.
% 5) Computes population-level summary statistics (mean and propagated error)
%    for Delta, k, and t0 parameters.
% 6) Generates combined box plots comparing Tr and Tb across stimulation
%    frequencies for:
%       - Response amplitude (Delta)
%       - Decay rate (k)
%       - Normalized half-decay time (t0 / 200)
% =========================================================================

% Directory containing .mat files with sigmoid fit results
dataDir = 'C:\Users\matte\Documents\Università\Magistrale\Secondo anno\Secondo semestre\TESI LM\Hyperione\NEW single electrodes';
matFiles = dir(fullfile(dataDir, "*Delta.mat"));

% Directory for saving output figures
figDir = fullfile(dataDir, 'figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

% Frequency tags as they appear in dataset names (human-readable)
freqTags = {'0_1Hz','0_2Hz','0_5Hz','1_0Hz','2_0Hz'};

% Valid MATLAB field names corresponding to the frequencies
validTags = {'Hz_0_1','Hz_0_2','Hz_0_5','Hz_1_0','Hz_2_0'};

% Structures used to collect filtered Tr and Tb results per frequency
resultsTr = struct();
resultsTb = struct();

for f = 1:numel(validTags)
    resultsTr.(validTags{f}) = table();
    resultsTb.(validTags{f}) = table();
end

%% ========================================================================
%  Aggregate fit results across all .mat files, grouping by stimulation frequency

% Loop over all .mat files in the directory
for fileIdx = 1:numel(matFiles)

    % Load the current file (must contain a variable named "fitResults")
    S = load(fullfile(matFiles(fileIdx).folder, matFiles(fileIdx).name));
    fitResults = S.fitResults;

    % List all dataset names stored inside this fitResults struct
    datasetNames = fieldnames(fitResults);

    % Loop over the target stimulation frequencies to be extracted
    for f = 1:numel(freqTags)
        tag = freqTags{f}; % string tag used in dataset names (e.g., '0_5Hz')
        fieldName = validTags{f}; % valid MATLAB field name used for storage (e.g., 'Hz_0_5')

        % Identify which datasets correspond to the current frequency tag
        matchIdx = contains(datasetNames, tag);
        matchedNames = datasetNames(matchIdx);

        % Loop over the datasets matching the current frequency
        for m = 1:numel(matchedNames)
            ds = fitResults.(matchedNames{m}); % structure containing Tr and Tb fit arrays

            % ---------------------------
            % Tr fits: append QC-passed rows
            for i = 1:numel(ds.Tr)
                row = ds.Tr(i);

                % QC criteria:
                % 1) Delta must be > 1e-3 (avoid near-flat/degenerate fits)
                % 2) maxBeforeMin == 1 (expected monotonic sigmoid direction)
                if all(row.Delta > 1e-3) && all(row.maxBeforeMin == 1)
                    resultsTr.(fieldName) = [resultsTr.(fieldName); struct2table(row)];
                end
            end

            % ---------------------------
            % Tb fits: append QC-passed rows
            for i = 1:numel(ds.Tb)
                row = ds.Tb(i);

                % Apply the same QC criteria used for Tr
                if all(row.Delta > 1e-3) && all(row.maxBeforeMin == 1)
                    resultsTb.(fieldName) = [resultsTb.(fieldName); struct2table(row)];
                end
            end

        end
    end
end

%% ========================================================================
% Remove duplicate rows across aggregated tables
for f = 1:numel(validTags)
    if ~isempty(resultsTr.(validTags{f}))
        resultsTr.(validTags{f}) = unique(resultsTr.(validTags{f}), 'rows');
    end
    if ~isempty(resultsTb.(validTags{f}))
        resultsTb.(validTags{f}) = unique(resultsTb.(validTags{f}), 'rows');
    end
end

%% ========================================================================
%  Manual removal of specific outlier rows identified during inspection
if height(resultsTr.Hz_1_0) >= 39
    resultsTr.Hz_1_0(39,:) = [];
end
if height(resultsTr.Hz_0_2) >= 3
    resultsTr.Hz_0_2(3,:) = [];
end
if height(resultsTr.Hz_0_1) >= 8
    resultsTr.Hz_0_1([7 8],:) = [];
elseif height(resultsTr.Hz_0_1) >= 7
    resultsTr.Hz_0_1(7,:) = [];
end

%% ========================================================================
%  ===== Summary statistics (mean + propagated error) =====

summaryTr = table();
summaryTb = table();

for f = 1:numel(validTags)
    tag = validTags{f};

    % --- Tr ---
    T = resultsTr.(tag);
    if ~isempty(T)
        N = height(T);

        meanDelta = mean(T.Delta);
        errDelta  = sqrt(sum(T.Delta_std.^2)) / N;

        kErrs = (T.k_ci_high_sec - T.k_ci_low_sec) / 2;
        meanK = mean(T.k_sec);
        errK  = sqrt(sum(kErrs.^2)) / N;

        t0Errs = (T.t0_ci_high - T.t0_ci_low) / 2;
        meanT0 = mean(T.t0 / 200);
        errT0  = sqrt(sum((t0Errs / 200).^2)) / N;

        summaryTr = [summaryTr;
            table({tag}, meanDelta, errDelta, meanK, errK, meanT0, errT0, ...
            'VariableNames', {'Freq','Delta','Delta_err','k','k_err','t0','t0_err'})];
    end

    % --- Tb ---
    T = resultsTb.(tag);
    if ~isempty(T)
        N = height(T);

        meanDelta = mean(T.Delta);
        errDelta  = sqrt(sum(T.Delta_std.^2)) / N;

        kErrs = (T.k_ci_high_sec - T.k_ci_low_sec) / 2;
        meanK = mean(T.k_sec);
        errK  = sqrt(sum(kErrs.^2)) / N;

        t0Errs = (T.t0_ci_high - T.t0_ci_low) / 2;
        meanT0 = mean(T.t0 / 200);                  
        errT0  = sqrt(sum((t0Errs / 200).^2)) / N;

        summaryTb = [summaryTb;
            table({tag}, meanDelta, errDelta, meanK, errK, meanT0, errT0, ...
            'VariableNames', {'Freq','Delta','Delta_err','k','k_err','t0','t0_err'})];
    end
end

%% ========================================================================
% labels for frequencies
freqLabels = {'0.1','0.2','0.5','1.0','2.0'};
nFreq = numel(freqLabels);
xShift = 0.15; 

%% ========================================================================
%  ===== Prepare data for boxplots =====

allDeltaTr = []; groupDeltaTr = [];
allKTr     = []; groupKTr     = [];
allT0Tr    = []; groupT0Tr    = [];

allDeltaTb = []; groupDeltaTb = [];
allKTb     = []; groupKTb     = [];
allT0Tb    = []; groupT0Tb    = [];

for f = 1:numel(validTags)
    if ~isempty(resultsTr.(validTags{f}))
        T = resultsTr.(validTags{f});
        g = f * ones(height(T),1);

        allDeltaTr = [allDeltaTr; T.Delta];
        groupDeltaTr = [groupDeltaTr; g];

        allKTr = [allKTr; T.k_sec];
        groupKTr = [groupKTr; g];

        allT0Tr = [allT0Tr; T.t0/200];
        groupT0Tr = [groupT0Tr; g];
    end

    if ~isempty(resultsTb.(validTags{f}))
        T = resultsTb.(validTags{f});
        g = f * ones(height(T),1);

        allDeltaTb = [allDeltaTb; T.Delta];
        groupDeltaTb = [groupDeltaTb; g];

        allKTb = [allKTb; T.k_sec];
        groupKTb = [groupKTb; g];

        allT0Tb = [allT0Tb; T.t0/200];
        groupT0Tb = [groupT0Tb; g];
    end
end

%% ========================================================================
%  ===== Combined box plots (Tr + Tb) =====

% Common style settings
colTr   = [0.6 0 0];
colTb   = [0 0 0.6];
symTr   = 'o';
symTb   = 's';
boxW    = 0.25;

% ---- Delta: plot only 0.5, 1.0, 2.0 Hz (groups 3..nFreq) but keep full x-axis
keepGroupsDelta = 3:nFreq;
plotBoxTrTb(allDeltaTr, groupDeltaTr, allDeltaTb, groupDeltaTb, ...
    nFreq, freqLabels, xShift, keepGroupsDelta, ...
    'Stimulation frequency [Hz]', '\Delta [1/s]', 'Firing rate decline (T_r vs T_b)', ...
    [0 24], [], {'T_r data','T_b data'}, colTr, colTb, symTr, symTb, boxW, ...
    fullfile(figDir, 'Delta_box_combined.png'));

% ---- k: plot all frequencies
plotBoxTrTb(allKTr, groupKTr, allKTb, groupKTb, ...
    nFreq, freqLabels, xShift, 1:nFreq, ...
    'Stimulation frequency [Hz]', 'k [1/s]', 'Decay rate (Tr vs Tb)', ...
    [0 0.16], [0 0.025 0.05 0.075 0.1 0.125 0.15], {'Tr','Tb'}, colTr, colTb, symTr, symTb, boxW, ...
    fullfile(figDir, 'k_box_combined.png'));

% ---- t0/200: plot all frequencies
plotBoxTrTb(allT0Tr, groupT0Tr, allT0Tb, groupT0Tb, ...
    nFreq, freqLabels, xShift, 1:nFreq, ...
    'Stimulation frequency [Hz]', 't_{0}/200', 'Fractional half-decay time (Tr vs Tb)', ...
    [], [], {'Tr','Tb'}, colTr, colTb, symTr, symTb, boxW, ...
    fullfile(figDir, 't0_box_combined.png'));

%% ========================================================================
%  ===== Local function for boxplots ======

function plotBoxTrTb(yTr, gTr, yTb, gTb, nFreq, freqLabels, xShift, keepGroups, ...
    xLab, yLab, ttl, yLims, yTicks, legTxt, colTr, colTb, symTr, symTb, boxW, outFile)
    % Plot combined Tr/Tb boxplots with optional group filtering 
    % (e.g., exclude low freqs) while keeping a full x-axis (1..nFreq)

    fig = figure; hold on;

    % Filter groups to be displayed (keeps axis complete regardless)
    idxTr = ismember(gTr, keepGroups);
    idxTb = ismember(gTb, keepGroups);

    % Positions must correspond to group numbers (so empty groups are left blank)
    posTr = keepGroups - xShift;
    posTb = keepGroups + xShift;

    % Tr boxplot
    boxplot(yTr(idxTr), gTr(idxTr), ...
        'Positions', posTr, 'Colors', colTr, 'Symbol', symTr, 'Widths', boxW);

    % Tb boxplot
    boxplot(yTb(idxTb), gTb(idxTb), ...
        'Positions', posTb, 'Colors', colTb, 'Symbol', symTb, 'Widths', boxW);

    % Aesthetic: thicker boxes/medians
    set(findobj(gca,'Tag','Box'),   'LineWidth', 1.2);
    set(findobj(gca,'Tag','Median'),'LineWidth', 1.2);

    % Full x-axis labels (even if some groups are excluded)
    set(gca, 'XTick', 1:nFreq, 'XTickLabel', freqLabels);
    xlabel(xLab);
    ylabel(yLab);
    title(ttl);
    grid on;
    axis square;

    % Optional y formatting
    if ~isempty(yLims),  ylim(yLims);  end
    if ~isempty(yTicks), yticks(yTicks); end

    % Legend (dummy handles)
    h1 = plot(nan, nan, 's', 'MarkerFaceColor', colTr, 'Color', colTr);
    h2 = plot(nan, nan, 's', 'MarkerFaceColor', colTb, 'Color', colTb);
    legend([h1 h2], legTxt, 'Location','best');

    saveas(fig, outFile);
    close(fig);
end
