% =========================================================================
% This script loads and aggregates ΔA (efficacy) measurements from 
% multiple .mat files corresponding to independent LTP experiments. 
% For each experiment, ΔA values from blue and red electrodes are collected 
% across multiple time points, including a pre-stimulation baseline.
%
% The script ensures consistency in data structure by adding a missing
% pre-stimulation column when necessary, then computes the global mean and
% variability (standard deviation or standard error of the mean) across all
% experiments.
%
% Finally, the script plots the average ΔA time course for blue and red
% electrodes, including error bars, to visualize the temporal evolution of
% synaptic changes following stimulation.
% =========================================================================

%% ========================================================================
%  ===== Data aquisition =====

% Directory containing the .mat files
dataDir = 'C:\Users\matte\Documents\Università\Pubblicazione Articolo\Codici\LTP_experiments';

% Get list of all .mat files in the directory
files = dir(fullfile(dataDir,'*.mat'));
if isempty(files)
    error('No .mat files found in %s', dataDir);
end

% Global accumulators for all experiments
Delta_blue_all = [];
Delta_red_all  = [];

%% ========================================================================
%  ===== Load and merge data from all files =====

for f = 1:numel(files)
    fname = fullfile(files(f).folder, files(f).name);
    fprintf('Loading %s\n', fname);
    D = load(fname);
    
    % Check and standardize matrix size
    % If the "pre" column is missing (only 6 columns), prepend a NaN column
    if isfield(D,'Delta_blue_all')
        Db = D.Delta_blue_all;
        if size(Db,2) == 6
            Db = [nan(size(Db,1),1), Db]; % add pre column
        end
        Delta_blue_all = [Delta_blue_all; Db];
    else
        warning('File %s does not contain Delta_blue_all', fname);
    end
    
    if isfield(D,'Delta_red_all')
        Dr = D.Delta_red_all;
        if size(Dr,2) == 6
            Dr = [nan(size(Dr,1),1), Dr]; % add pre column
        end
        Delta_red_all = [Delta_red_all; Dr];
    else
        warning('File %s does not contain Delta_red_all', fname);
    end
end

%% ========================================================================
%  ===== Definition of Time Points (min) =====

% Includes the pre-stimulation baseline
post_times = [-37 0 30 60 90 120 180];

%% ========================================================================
%  ===== Mean and variability estimation =====

% Compute mean, standard deviation, and standard error of the mean (SEM)
% both for blue and red electrodes
Delta_blue = mean(Delta_blue_all, 1, 'omitnan');
std_blue   = std(Delta_blue_all, 0, 1, 'omitnan');
n_blue     = size(Delta_blue_all, 1);
sem_blue   = std_blue ./ sqrt(n_blue);

Delta_red = mean(Delta_red_all, 1, 'omitnan');
std_red   = std(Delta_red_all, 0, 1, 'omitnan');
n_red     = size(Delta_red_all, 1);
sem_red   = std_red ./ sqrt(n_red);

%% ========================================================================
%  ===== PLOT INCLUDING ALL TIME POINTS (PRE INCLUDED) =====

figure; hold on;

% Choose whether to plot SEM or STD as error bars
use_sem = true;
if use_sem
    err_blue = sem_blue;
    err_red  = sem_red;
    err_label = 'SEM';
else
    err_blue = std_blue;
    err_red  = std_red;
    err_label = 'STD';
end

% Dashed line with markers for Blue electrodes
errorbar(post_times, Delta_blue, err_blue, '--o', ...
    'Color',[0 0 1], 'LineWidth',1.5, 'CapSize',5, ...
    'MarkerFaceColor',[0.3 0.3 1]);

% Dashed line with markers for Red electrodes
errorbar(post_times, Delta_red, err_red, '--o', ...
    'Color',[1 0 0], 'LineWidth',1.5, 'CapSize',5, ...
    'MarkerFaceColor',[1 0.3 0.3]);

xlabel('Time (min)');
ylabel('Average \DeltaA');
legend('Blue electrodes','Red electrodes', 'Location','best');
title('Global average \DeltaA as a function of time');
grid on;

% X-axis formatting
xlim([min(post_times) max(post_times)]);
xticks(post_times);
xticklabels(arrayfun(@num2str, post_times, 'UniformOutput', false));

% Automatic y-axis limits including error bars
Ymin = min([Delta_blue - err_blue, Delta_red - err_red], [], 'omitnan');
Ymax = max([Delta_blue + err_blue, Delta_red + err_red], [], 'omitnan');
if isnan(Ymin) || isnan(Ymax)
    ylim([0 1]);
elseif Ymax == Ymin
    ylim([Ymin - 0.05, Ymax + 0.05]);
else
    rangeY = Ymax - Ymin;
    ylim([Ymin - 0.1*rangeY, Ymax + 0.1*rangeY]);
end

% Improve readability
set(gca,'FontSize',11);
