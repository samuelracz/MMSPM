function [ gfits ] = func_Gaussian_peaks(Y, freq, m_set)

f_options = optimoptions('fmincon','Display','off');

Y_tmp = Y;
th_Amp = m_set.min_peak_Amp;
th_std = m_set.min_peak_SD;
max_N = m_set.max_peak_N;
gfits = [];

while true
    % if maximum number of peaks surpassed, break
    if length(gfits)/3 > max_N
        break
    end

    % spectrum parameters
    Y_mu = mean(Y_tmp);
    Y_std = std(Y_tmp);

    % detect peak
    [A_max, f_max] = max(Y_tmp);
    g_mu = freq(f_max);

    % find FWHM+ and FWHM- around peak
    [~,fwhm_n,fwhm_p] = func_fwhm(freq,Y_tmp);
    fwhm = min([fwhm_n, fwhm_p]);
    if fwhm == 0
        fwhm = min(diff(freq));
    end
    g_std = fwhm/(2*sqrt(2*log(2)));

    % constrained optimization for Gaussian curve fitting
    cost_Gauss = @(x)sum((Y_tmp - (x(1)*exp(-(1/2)*(((freq-x(2))/(x(3))).^2)))).^2);
    x0 = [A_max; g_mu; g_std];
    lb = [0; g_mu-2*fwhm; 0];
    ub = [Y_mu+100*Y_std; g_mu+2*fwhm; 2*fwhm];

    gf = fmincon(cost_Gauss,x0,[],[],[],[],lb,ub,[],f_options);
    gcurve = gf(1)*exp(-(1/2)*(((freq-gf(2))/(gf(3))).^2));

    % plot for sanity check
    if m_set.debug
        f = figure;
        plot(freq,Y_tmp)
        hold on
        plot([min(freq) max(freq)], [Y_mu + th_std*Y_std Y_mu + th_std*Y_std],'k--')
        plot([min(freq) max(freq)], [th_Amp th_Amp],'k.-')
        plot(freq,gcurve,':','LineWidth',2)
        plot([freq(f_max) freq(f_max)],[0, A_max],'go-','LineWidth',2)
        plot([freq(f_max) freq(f_max)-fwhm_n],[A_max/2, A_max/2],'gd-','LineWidth',2)
        plot([freq(f_max) freq(f_max)+fwhm_p],[A_max/2, A_max/2],'gd-','LineWidth',2)
        input('press ENTER to continue')
        close(f)
    end

    % if no peak, break loop
    if gf(1) < Y_mu + th_std*Y_std
        break
    end

    if gf(1) < th_Amp
        break
    end

    % store gaussian fit
    if isempty(gfits)
        gfits = gf(:);
    else
        gfits_tmp = cat(1,gfits,gf(:));
        gfits = gfits_tmp;
    end

    % remove gaussian peak
    Y_tmp = Y_tmp - gcurve;

end

end

%% helper functions
function [fm, fwhm_neg, fwhm_pos] = func_fwhm(freq,Y)
% find full widht at half maximum in both directions

Nf = length(freq);

% maximum power
[fm, max_ind] = max(Y);

% FWHM-
fwhm_neg = [];
ind = max_ind;
while ind > 1
    ind = ind-1;
    % first time power drops below fm/2
    if Y(ind) < fm/2
        fwhm_neg = freq(max_ind) - freq(ind);
        break
    end
end

% if empty
if isempty(fwhm_neg)
    [~, min_ind] = min(Y(1:max_ind));
    fwhm_neg = freq(max_ind) - freq(min_ind);
end

% FWHM+
fwhm_pos = [];
ind = max_ind;
while ind < Nf
    ind = ind+1;
    % first time power drops below fm/2
    if Y(ind) < fm/2
        fwhm_pos = freq(ind) - freq(max_ind);
        break
    end
end

% if empty
if isempty(fwhm_pos)
    [~, min_ind] = min(Y(max_ind:end));
    fwhm_pos = freq(max_ind+min_ind-1) - freq(max_ind);
end

end

