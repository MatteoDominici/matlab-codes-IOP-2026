function [t, Data, Channels] = AnalogDataExtr(FileName, n1, n2)
    % =====================================================================
    % The function extracts analogue signals, time vector, and channel 
    % identifiers from an HDF5 recording file.
    %
    % INPUTS:
    %   - FileName (string) : path to the HDF5 recording file
    %   - n1 (integer) : index of the first sample to extract
    %   - n2 (integer) : number of samples to extract
    %
    % OUTPUTS:
    %   - t (vector) : time vector corresponding to the extracted samples (s)
    %   - Data ([nSamples x nChannels] matrix) : analogue signal data expressed in µV
    %   - Channels (vector) : numeric channel identifiers
    % =====================================================================
    
    % ===== Channel metadata extraction =====
    % Read channel information from the analogue stream
    Info = h5read(FileName, ...
        '/Data/Recording_0/AnalogStream/Stream_0/InfoChannel');

    % Determine the number of channels
    Nch = length(Info.Tick);

    % ===== Analogue data extraction =====
    % Read analogue signal data for the selected sample range and channels
    Data = double(h5read(FileName, ...
        '/Data/Recording_0/AnalogStream/Stream_0/ChannelData', [n1 1], [n2 Nch]));

    % ===== Unit conversion to microvolts =====
    % Retrieve conversion parameters
    Factor = double(Info.ConversionFactor(1));
    EXP    = 10 ^ double(Info.Exponent(1));

    % Convert raw data to microvolts (µV)
    Data = 1e6 * Data * Factor * EXP;

    % ===== Time vector construction =====
    % Generate sample index vector
    cnt = (0:size(Data,1)-1)';

    % Convert sample indices to time (s)
    t = 1e-6 * double(Info.Tick(1)) * cnt;

    % ===== Channel label extraction =====
    % Convert channel labels to numeric format
    Channels = str2double(Info.Label);
end
