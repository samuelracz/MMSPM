function [ Y ] = func_generate_power_spectrum(freq, params)

N = length(freq);

b_lo = -1*params.b_lo;      % low-range slope
b_hi = -1*params.b_hi;      % high-range slope
bp = params.bp;             % breakpoint (in Hz)
trw = params.trw;           % transition range width (in Hz)
snr = params.snr;           % noise level (signal-to-noise ratio)
peaks = params.peaks;       % oscillatory peaks (center, spread, amplitude)

Np = size(peaks,1);

ind_lo = freq <= bp;
ind_hi = freq > bp;
xb = max(freq(ind_lo));


% obtain initial bimodal model
a1 = 10;
a2 = a1*(xb^b_lo)/(xb^b_hi);
Y_raw = log((a1*freq.^(b_lo))).*ind_lo + log((a2*freq.^(b_hi))).*ind_hi;

% add transition zone
if trw ~= 0
    ind_trw_lo = freq <= xb - 0.5*trw;
    ind_trw_hi = freq >= xb + 0.5*trw;
    ind_trw_out = boolean(ind_trw_lo + ind_trw_hi);
    ind_trw_in = boolean(-1*(ind_trw_lo + ind_trw_hi)+1);

    freq_trw_out = freq(ind_trw_out);
    freq_trw_in = freq(ind_trw_in);
    Y_trw_out = Y_raw(ind_trw_out);
    Y_trw_in = spline(log(freq_trw_out), Y_trw_out, log(freq_trw_in));

    Y_raw(ind_trw_in) = Y_trw_in;
end

Y_raw = Y_raw(:);
Y_raw = exp(Y_raw)./sum(exp(Y_raw));
Y_raw = log(Y_raw);

% add oscillatory peaks
Y_peaks = zeros(size(Y_raw));
if ~isempty(peaks)
    for n = 1:Np
        cf = peaks(n,1);
        ps = peaks(n,2);
        pA = peaks(n,3);

        Y_peaks = Y_peaks + pA*pdf(gmdistribution(cf,ps^2),freq)./pdf(gmdistribution(cf,ps^2),cf);
        % Y_peaks = Y_peaks + pA*pdf(gmdistribution(cf,ps^2),freq);
    end
end


% add noise
Y_noise = randn(N,1);
Y_noise = exp(Y_noise)./sum(exp(Y_noise));
Y_noise = log(Y_noise);

% full model
Y = exp(Y_raw + Y_peaks + snr.*Y_noise);




end