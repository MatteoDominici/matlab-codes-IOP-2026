function [psth_norm, area_spikes] = ComputePSTH(spikes, tpOn, binSize, win)
    % =====================================================================
    % The function computes the peri-stimulus time histogram (PSTH) aligned 
    % to stimulation onset times.
    %
    % INPUTS:
    %   - spikes (vector) : spike times (s)
    %   - tpOn (vector) : stimulation onset times (s)
    %   - binSize (double) : width of each PSTH bin (s)
    %   - win (1x2 vector) : time window relative to stimulation onset (s)
    %
    % OUTPUTS:
    %   - psth_norm (vector) : normalized PSTH expressed as firing rate (Hz)
    %   - area_spikes (double) : total spike count within the window, averaged per stimulus
    % =====================================================================
    
    % ===== PSTH bin definition =====
    % Define bin edges spanning the analysis window
    edges = win(1):binSize:win(2);

    % Determine the number of histogram bins
    nBins = length(edges) - 1;

    % Initialize PSTH count accumulator
    psth_counts = zeros(1, nBins);

    % ===== Spike alignment and Binning =====
    % Loop over all stimulation onset times
    for i = 1:length(tpOn)

        % Current stimulation onset time
        stim_time = tpOn(i);

        % Absolute start and end of the analysis window
        win_start = stim_time + win(1);
        win_end   = stim_time + win(2);

        % Extract spikes occurring within the analysis window
        % and align them to stimulation onset (time zero)
        spikes_in_window = spikes(spikes >= win_start & ...
                                  spikes <  win_end) - stim_time;

        % Skip iteration if no spikes are present in the window
        if isempty(spikes_in_window), continue; end

        % Compute histogram counts for the current stimulus
        counts = histcounts(spikes_in_window, edges);

        % Accumulate counts across stimuli
        psth_counts = psth_counts + counts;
    end

    % ===== PSTH normalization and spike area computation =====
    % Normalize PSTH by number of stimuli and bin width (firing rate in Hz)
    psth_norm = psth_counts ./ length(tpOn) ./ binSize;

    % Compute the total spike area under the PSTH curve
    area_spikes = sum(psth_norm) * binSize;
end
