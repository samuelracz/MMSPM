function [ settings_out ] = func_set_mmspm_settings(settings_in)

%% set defaults
settings_out = struct(...
     's_env', [1, 2, 4, 8, 16],...          % window sizes for envelope fitting (Hz)
     'f_bp', 'auto',...                     % initial guess for breakpoint: either 'auto' or [initial guess, lower limit, upper limit]
     'f_res', [],...                        % frequency resolution - leave empty
     'f_ma', 1,...                          % moving average smoothing window size (Hz)
     'f_w', 1,...                           % exponent for weighting frequencies to compensate for uneven sampling
     'q_lo', 0.25,...                       % lower quantile for initial bimodal fit
     'max_peak_N', 10,...                   % maximum number of peaks
     'min_peak_Amp', 0,...                  % minimum peak amplitude in power
     'min_peak_SD', 1.2,...                 % minimum peak amplitude in terms of spectrum standard deviation
     'min_peak_detect_width', 1,...         % minimum peak width in peak detection
     'th_peak_width', 2,...               % threshold for excluding peak neighborhood (initial model guess)
     'th_peak_height', 0,...                % threshold for excluding peak neighborhood (initial model guess)
     'th_edge', 1,...                       % threshold for excluding peaks close to analysis range (in SD)
     'th_overlap', 0.25,...                 % threshold for excluding overlapping peaks (in SD)
     'th_dist', 1.5,...                     % threshold distance from original peak guess in Gaussian fitting
     'th_exclude', 0,...                    % threshold in excluding frequency ranges with peaks from bimodal fit
     'bimodal_alpha', 0.05,...              % statistical significance threshold in testing for true bimodality
     'flag_fig', true,...                   % flag for plotting outcomes
     'flag_disp', true,...                  % flag for displaying results in command window
     'debug', false);                       % debug mode

%% overwrite defaults based on input
vnames = fieldnames(settings_out);
vnames_in = fieldnames(settings_in);

for v = 1:length(vnames_in)
    if contains(vnames_in{v},vnames)
        settings_out.(vnames_in{v}) = settings_in.(vnames_in{v});
    end
end

end