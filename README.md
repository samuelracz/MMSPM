# Multi-Modal Spectral Parametrization Method (MMSPM) for analyzing EEG activity with distinct scaling regimes

## This software is still in development.

The purpose of MMSPM is to characterize the power spectra of brain electrophysiological recordings (such as EEG or ECoG) by deconstructing it into broadband aperiodic and narrow-band oscillatory components. In this MMSPM framework, the broadband component is modeled as a composite of two power-law fucntions with respective offsets _A<sub>1</sub>_ and _A<sub>2</sub>_ and scaling exponents _b<sub>1</sub>_ and _b<sub>2</sub>_; the two power-law functions have independent scaling regimes separated by a breakpoint frequency _f<sub>bp</sub>_, and the constraint _A<sub>1</sub>*f<sub>bp</sub><sup>b<sub>1</sub></sup>=A<sub>2</sub>*f<sub>bp</sub><sup>b<sub>2</sub></sup>_ is enforced assuming continuity of the aperiodic component. Superimposed oscillatory peaks are modeled as independent Gaussians.

## Usage:
First, a parameter settings Matlab structure needs to be defined. More details on default parameter settings can be found in _func_set_mmspm_settings.m_. The list of possible parameter settings are:
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
    'th_peak_width', 2,...                 % threshold for excluding peak neighborhood (initial model guess)
    'th_peak_height', 0,...                % threshold for excluding peak neighborhood (initial model guess)
    'th_edge', 1,...                       % threshold for excluding peaks close to analysis range (in SD)
    'th_overlap', 0.25,...                 % threshold for excluding overlapping peaks (in SD)
    'th_dist', 1.5,...                     % threshold distance from original peak guess in Gaussian fitting
    'th_exclude', 0,...                    % threshold in excluding frequency ranges with peaks from bimodal fit
    'bimodal_alpha', 0.05,...              % statistical significance threshold in testing for true bimodality
    'flag_fig', true,...                   % flag for plotting outcomes
    'flag_disp', true,...                  % flag for displaying results in command window
    'debug', false,...                     % debug mode

    
