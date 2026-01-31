function saveSigmoidFitPlot(x, y_raw, y_smooth, y_fit, fname, chID, R2, tag, outDir)
% =========================================================================
% Generates and saves a diagnostic plot for a sigmoid fit of spike-rate data.
% The figure displays:
%   - Raw spike-rate values (per pulse)
%   - Smoothed spike-rate time course
%   - Fitted sigmoid model
% The plot is exported as a high-resolution .jpg file and then closed to
% avoid accumulation of open figures during batch processing.
%
% INPUTS:
%   - x        : pulse indices (independent variable)
%   - y_raw    : raw spike-rate values
%   - y_smooth : smoothed spike-rate values used for fitting
%   - y_fit    : fitted sigmoid evaluated at x
%   - fname    : base filename (used in title and output filename)
%   - chID     : MEA channel ID
%   - R2       : coefficient of determination of the fit
%   - tag      : string label identifying condition ('Tr', 'Tb', etc.)
%   - outDir   : output directory for saved figures
%
% OUTPUT:
%   - None (figure is saved to disk)
% =========================================================================

    % Create figure with standardized size and white background
    fig = figure( ...
        'Units','normalized', ...
        'OuterPosition',[0 0 0.5 0.5], ...
        'Color','w');

    % Plot raw data (markers only)
    plot(x, y_raw, 'x', ...
        'Color',[0.25 0.25 0.25], ...
        'MarkerSize', 6); 
    hold on;

    % Plot smoothed spike-rate trace
    plot(x, y_smooth, 'b-', 'LineWidth', 1.2);

    % Plot fitted sigmoid curve
    plot(x, y_fit, 'r-', 'LineWidth', 1.5);

    % Axis labels
    xlabel('Pulse #', 'FontSize', 10);
    ylabel('Rate (Hz)', 'FontSize', 10);

    % Title with metadata and goodness-of-fit
    title(sprintf('Sigmoid fit %s — %s — Ch %d (R^2=%.2f)', ...
        tag, fname, chID, R2), ...
        'Interpreter','none', ...
        'FontSize', 11);

    % Legend
    legend({'Raw','Smoothed','Fit'}, ...
        'FontSize', 9, ...
        'Location','best');

    % Axis formatting (fixed limits for comparability across channels)
    ylim([0 30]);
    yticks(0:5:30);
    axis square;
    set(gca,'FontSize', 9, 'Box','on');
    grid on;

    % Create output directory if it does not exist
    if ~exist(outDir,'dir')
        mkdir(outDir);
    end

    % Construct output filename and export figure
    outName = sprintf('%s_%s_fit_ch%d.jpg', fname, tag, chID);
    exportgraphics(fig, fullfile(outDir, outName), 'Resolution', 300);

    % Close figure to prevent memory/handle accumulation
    close(fig);
end
