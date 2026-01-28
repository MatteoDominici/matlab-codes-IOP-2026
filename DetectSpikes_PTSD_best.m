function spike_times = DetectSpikes_PTSD_best(filtered, fs, k)
    % =====================================================================
    % The function detects spike events from band-pass filtered 
    % electrophysiological signals using a peak-to-peak thresholding 
    % strategy with temporal constraints.
    %
    % INPUTS:
    %   - filtered ([nSamples x nChannels] matrix) : filtered signal
    %   - fs (double) : sampling frequency (Hz)
    %   - k (double) : threshold multiplier (DT = k * sigma)
    %
    % OUTPUT:
    %   - spike_times (1 x nChannels cell array) : each cell contains spike 
    %               sample indices for the corresponding channel
    % =====================================================================

    % ===== Temporal parameters (in milliseconds) ======
    PLP_ms       = 1.5;  % Peak Lifetime Period
    RP_ms        = 1.0;  % Refractory Period
    overshoot_ms = 0.5;  % Additional search window if peak lies at PLP edge

    % ===== Conversion from milliseconds to samples =====
    PLP       = round(PLP_ms       * 1e-3 * fs);
    RP        = round(RP_ms        * 1e-3 * fs);
    overshoot = round(overshoot_ms * 1e-3 * fs);

    % ===== Initialization =====
    [nSamples, nChannels] = size(filtered);
    spike_times = cell(1, nChannels);

    % ===== Channel-wise spike detection =====
    for ch = 1:nChannels

        % Extract signal for the current channel
        sig = filtered(:, ch);

        % Skip channels containing only zeros or NaN values
        if all(sig == 0) || any(isnan(sig))
            fprintf('Skipping channel %d: flat signal or NaN values detected\n', ch);
            continue;
        end

        % ===== Baseline statistics estimation =====
        % Estimate mean and standard deviation from the first second of data
        baseN = min(fs, nSamples);
        mu    = mean(sig(1:baseN));
        sigma = std(sig(1:baseN));
        DT = k * sigma; % Differential threshold

        % ===== Relative minima and maxima detection =====
        isMin = islocalmin(sig, 1);
        isMax = islocalmax(sig, 1);

        % Identify noisy extrema:
        % - Low maxima: maxima below the mean
        % - High minima: minima above the mean
        low_maxima  = false(nSamples,1);
        high_minima = false(nSamples,1);
        if any(isMax), low_maxima(isMax)  = sig(isMax) < mu; end
        if any(isMin), high_minima(isMin) = sig(isMin) > mu; end

        % ===== Candidate extrema selection =====
        % Retain only extrema compatible with spike morphology
        candMask = (isMin & ~high_minima) | (isMax & ~low_maxima);
        candIdx  = find(candMask);

        % Skip channel if no candidates are found
        if isempty(candIdx)
            fprintf('Channel %d: no candidate extrema found\n', ch);
            continue;
        end

        % Pointer array used to mark samples already occupied by
        % refractory periods or previously paired extrema
        pointer = false(nSamples,1);

        % Container for detected spike indices
        spikes = [];

        % ===== Candidate processing loop =====
        for t0 = candIdx.'   % t0 is a candidate relative minimum or maximum

            % Skip candidates too close to signal end or already blocked
            if t0 + PLP >= nSamples, break; end
            if pointer(t0), continue; end

            % Determine whether the candidate is a minimum or a maximum
            isT0Min = isMin(t0);

            % Define PLP search window
            winEnd = min(t0 + PLP, nSamples);

            % Search for the opposite extremum within PLP
            if isT0Min
                [~, pkRel] = max(sig(t0:winEnd));
            else
                [~, pkRel] = min(sig(t0:winEnd));
            end
            pkIdx = t0 + pkRel - 1;

            % Extend search window if the opposite peak lies on the PLP edge
            if pkIdx == winEnd && (winEnd < nSamples) && overshoot > 0
                newEnd = min(t0 + PLP + overshoot, nSamples);
                if isT0Min
                    [~, pkRel] = max(sig(t0:newEnd));
                else
                    [~, pkRel] = min(sig(t0:newEnd));
                end
                pkIdx = t0 + pkRel - 1;
            end

            % Skip if the opposite peak is already blocked
            if pointer(pkIdx), continue; end

            % ===== Peak-to-peak amplitude test =====
            segLo = min(t0, pkIdx);
            segHi = max(t0, pkIdx);
            seg   = sig(segLo:segHi);

            % Peak-to-peak amplitude
            p2p = max(seg) - min(seg);
            if p2p < DT, continue; end

            % ===== Spike timestamp selection =====
            % Collect all extrema within the segment
            segIdxAbs = (segLo:segHi).';
            segMins   = isMin(segIdxAbs);
            segMaxs   = isMax(segIdxAbs);
            rmmSeg    = segIdxAbs(segMins | segMaxs);

            % Ensure endpoints are included
            rmmSeg = unique([rmmSeg; t0; pkIdx]);

            % Select the extremum with maximum absolute deviation from the mean
            [~, bestPos] = max(abs(sig(rmmSeg) - mu));
            spikeIdx = rmmSeg(bestPos);

            % ===== Refractory period and segment blocking =====
            rpEnd = min(spikeIdx + RP, nSamples);
            pointer(spikeIdx:rpEnd) = true;
            pointer(segLo:segHi)    = true;
            
            spikes(end+1) = spikeIdx; % Store detected spike
        end

        % ===== Output assignment and reporting =====
        fprintf('Channel %d: %d spikes detected\n', ch, numel(spikes));
        spike_times{ch} = spikes;
    end
end
