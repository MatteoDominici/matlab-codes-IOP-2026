function fitStruct = fitSigmoidForSelectedChannels( ...
    Channels, SR_raw, SR_smooth, selectedChannels, ...
    T, outDir, fname, tag)
% =========================================================================
% Fits a 2-parameter sigmoid model to spike-rate time courses for selected
% MEA channels, computes goodness-of-fit metrics and confidence intervals,
% converts parameters to physical units (seconds), and saves diagnostic plots.
%
% INPUTS:
%   - Channels         : vector of MEA channel IDs
%   - SR_raw           : raw spike-rate matrix (channels x pulses)
%   - SR_smooth        : smoothed spike-rate matrix (channels x pulses)
%   - selectedChannels : vector of channel IDs to be fitted
%   - T                : stimulation period (s)
%   - outDir           : output directory for figures
%   - fname            : base filename (used in plot titles and filenames)
%   - tag              : string label ('Tr' or 'Tb')
%
% OUTPUT:
%   - fitStruct        : structure containing fit results for each channel
% =========================================================================

fitStruct = struct;

% Initial guess and fitting options
k0 = 1/200;
opts = optimoptions('lsqcurvefit','Display','off');
warning('off','stats:nlparci:JacobianNotFullRank');

for ch = 1:numel(Channels)

    thisCh = Channels(ch);

    % Fit only the selected MEA channels
    if ~ismember(thisCh, selectedChannels)
        continue;
    end

    y_smooth = SR_smooth(ch,:);
    y_raw    = SR_raw(ch,:);

    if isempty(y_smooth) || all(isnan(y_smooth))
        continue;
    end

    % Pulse index as independent variable
    y = y_smooth(:);
    x = (1:numel(y))';

    % Sigmoid asymptotes
    L = min(y);
    U = max(y);
    Delta = U - L;

    % Check ordering of extrema
    [~, idxMin] = min(y);
    [~, idxMax] = max(y);
    maxBeforeMin = (idxMax < idxMin);

    % Local uncertainty on Delta
    win = 2;
    rangeMin = max(1,idxMin-win):min(numel(y),idxMin+win);
    rangeMax = max(1,idxMax-win):min(numel(y),idxMax+win);
    stdL = std(y(rangeMin),'omitnan');
    stdU = std(y(rangeMax),'omitnan');
    stdDelta = sqrt(stdL^2 + stdU^2);

    % Sigmoid model
    t0_init = median(x);
    lb = [0 0];
    ub = [Inf numel(x)];
    sigmoid = @(b,xx) L + (U-L) ./ (1 + exp(b(1)*(xx-b(2))));

    try
        [beta,~,resid,~,~,~,J] = lsqcurvefit( ...
            sigmoid,[k0 t0_init],x,y,lb,ub,opts);

        k  = beta(1);
        t0 = beta(2);

        % Goodness of fit
        SSres = sum(resid.^2);
        SStot = sum((y-mean(y)).^2);
        R2 = 1 - SSres/SStot;

        % Confidence intervals
        ci = nlparci(beta,resid,'jacobian',J);

        % Conversion to seconds
        fitStruct(ch).channel        = thisCh;
        fitStruct(ch).k              = k;
        fitStruct(ch).k_ci_low       = ci(1,1);
        fitStruct(ch).k_ci_high      = ci(1,2);
        fitStruct(ch).t0             = t0;
        fitStruct(ch).t0_ci_low      = ci(2,1);
        fitStruct(ch).t0_ci_high     = ci(2,2);
        fitStruct(ch).R2             = R2;
        fitStruct(ch).Delta          = Delta;
        fitStruct(ch).Delta_std      = stdDelta;
        fitStruct(ch).maxBeforeMin   = maxBeforeMin;

        fitStruct(ch).t0_sec         = t0 * T;
        fitStruct(ch).t0_ci_low_sec  = ci(2,1) * T;
        fitStruct(ch).t0_ci_high_sec = ci(2,2) * T;
        fitStruct(ch).k_sec          = k / T;
        fitStruct(ch).k_ci_low_sec   = ci(1,1) / T;
        fitStruct(ch).k_ci_high_sec  = ci(1,2) / T;

        % Plot and save
        saveSigmoidFitPlot(x, y_raw(:), y, sigmoid(beta,x), ...
            fname, thisCh, R2, tag, outDir);

    catch ME
        warning('Fit failed for %s (%s, ch %d): %s', ...
            fname, tag, thisCh, ME.message);
    end
end
end
