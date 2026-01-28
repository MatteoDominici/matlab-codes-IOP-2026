function [spike_times_sec, Channels, tpOn] = AnalysisFile(filename, fs, k, bpf)
    % =====================================================================
    % The function performs spike detection and stimulation event extraction 
    % from an HDF5 recording file.
    %
    % INPUTS:
    %   - filename (string) : path to the HDF5 recording file
    %   - fs (double) : sampling frequency (Hz)
    %   - k (double) : threshold parameter for spike detection
    %   - bpf (1x2 vector) : band-pass filter cutoff frequencies (Hz)
    %
    % OUTPUTS:
    %   - spike_times_sec (cell array) : detected spike times (s) for each channel
    %   - Channels (array) : channel metadata extracted from the recording
    %   - tpOn (vector) : stimulation onset times (s) 
    %                     Empty if no stimulation events are present
    % =====================================================================

    fprintf('>>> Spike detection file: %s\n', filename);
    
    % ===== Analogue signal extraction =====
    % Retrieve information about the analogue data stream
    info = h5info(filename, '/Data/Recording_0/AnalogStream/Stream_0/ChannelData');

    % Extract analogue data and channel information
    [~, Data, Channels] = AnalogDataExtr(filename, 1, info.Dataspace.Size(1));

    % Ensure that data are organized as samples x channels
    if size(Data,1) < size(Data,2), Data = Data'; end

    % ===== Band-pass filtering =====
    % Design a third-order Butterworth band-pass filter
    [b,a] = butter(3, bpf/(fs/2), 'bandpass');

    % Apply zero-phase digital filtering to avoid phase distortion
    filtered = filtfilt(b,a,double(Data));

    % ===== Spike detection =====
    % Detect spike sample indices for each channel
    spike_samples = DetectSpikes_PTSD_best(filtered, fs, k);

    % Convert spike sample indices to time (seconds)
    spike_times_sec = cellfun(@(x) x/fs, spike_samples, 'UniformOutput', false);

    % ===== Stimulation event extraction (if available) =====
    % Initialize stimulation onset vector
    tpOn = [];

    try
        % Attempt to read stimulation ON events from the EventStream
        tpOn_raw = double(h5read(filename, ...
            '/Data/Recording_0/EventStream/Stream_0/EventEntity_0'));

        % Convert timestamps from microseconds to seconds
        tpOn = tpOn_raw(:,1) * 1e-6;

        % Attempt to read stimulation OFF events
        tpOff_raw = double(h5read(filename, ...
            '/Data/Recording_0/EventStream/Stream_1/EventEntity_0'));
        tpOff = tpOff_raw(:,1) * 1e-6;

        % Read analogue data timestamp offset (first sample)
        Dtp = 1e-6 * double(h5read(filename, ...
            '/Data/Recording_0/AnalogStream/Stream_0/ChannelDataTimeStamps'));
        Dtp = Dtp(1);

        % Align stimulation timestamps to analogue data start
        tpOn  = tpOn  - Dtp;
        tpOff = tpOff - Dtp;

        % Estimate stimulation period
        Tperiod = mode(diff(tpOn));

        % Remove last stimulation if it exceeds the recording length
        if (tpOn(end) + Tperiod) * fs > length(Data)
            tpOn(end)  = [];
            tpOff(end) = [];
        end

    catch ME
        % If no EventStream is present, assume spontaneous recording
        warning(['No stimulation events found in "%s". ' ...
                 'Assuming spontaneous activity.\nDetails: %s'], ...
                 filename, ME.message);

        % Explicitly return an empty stimulation onset vector
        tpOn = [];
    end
end
