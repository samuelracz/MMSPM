function [ o_power, o_peak, o_center ] = func_obtain_oscillatory_power(mmspm,frange)

freq = mmspm.freq;
peaks = mmspm.peaks;

if isempty(peaks)
    o_power = 0;
    o_peak = 0;
    o_center = 0;
else
    p_ind = peaks(:,2) > min(frange) & peaks(:,2) < max(frange);

    % select peaks of interest
    peaks = peaks(p_ind,:);
    Np = size(peaks,1);

    if Np ~= 0
        o_power_curve = zeros(length(freq),1);
        for n = 1:Np
            pA = peaks(n,1);
            cf = peaks(n,2);
            ps = peaks(n,3);

            o_power_curve = o_power_curve + pA*pdf(gmdistribution(cf,ps^2),freq)./pdf(gmdistribution(cf,ps^2),cf);
        end

        % integrate
        o_power = sum(o_power_curve);

        % select largest peak
        o_peak = max(o_power_curve);
        o_center = freq(o_power_curve==o_peak);
    else
        o_power = 0;
        o_peak = 0;
        o_center = 0;
    end
end

end