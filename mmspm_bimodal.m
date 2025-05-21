function [ b_struct ] = mmspm_bimodal(Y, freq, mmspm_settings)

% Multi-Modal Spectral Parametrization Method (MMSPM) for bimodal characterization of electrophysiological spectra
% Copyright (C) 2025 Frigyes Samuel Racz, Department of Neurology, Dell Medical School, The University of Texas at Austin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% DISCLOSURE: This program is still under development and therefore
% future updates to the code can be expected.
%
% The license applies equivalently to all functions called by mmspm_bimodal.m
%
% Last modified: 05/21/2025

%% Step 0. set defaults
m_set = func_set_mmspm_settings(mmspm_settings);

% frequency resolution
f_res = uniquetol(diff(freq),1e-5);
m_set.f_res = f_res;

% fmincon options
f_options = optimoptions('fmincon','Display','off');

%% Step 1. perform initial manipulations
% store for backup
Y_raw = Y;

% convert power to log-power
Y = log(Y);

% smooth with moving average
% Y = movavg(Y,'exponential',round(m_set.f_ma/f_res));

%% Step 1. obtain lower envelope for initial detrending
% W = round((1/f_res)*m_set.s_env);
% E_lo = min(func_env_low(Y, W),[],2);
E_lo = func_env_min(Y, freq, m_set);

%% Step 2. select points with lowest power
ind0 = func_env_nopeak(Y,freq,m_set.min_peak_detect_width,m_set.th_peak_height,m_set.th_peak_width);
Yd0 = Y - E_lo;

% add a small amount of noise to avoid zeros
Yd0 = Yd0 + rand(size(Yd0))*1e-10;

if isempty(ind0)
    ind0 = Yd0<=quantile(Yd0,m_set.q_lo);
end

%% Step 3. set 'initial guess' for breakpoint
logp = polyfit([log(freq(1)), log(freq(end))],[E_lo(1), E_lo(end)],1);
Y_logfit = log(exp(logp(2))*freq.^logp(1));
if strcmp(m_set.f_bp,'auto')
    [~, dind] = max(abs(E_lo - Y_logfit));
    bp_init = freq(dind);
    bp_lo = min(freq) + 0.05*range(freq);
    bp_hi = max(freq) - 0.05*range(freq);
else
    bp_init = m_set.f_bp(1);
    bp_lo = m_set.f_bp(2);
    bp_hi = m_set.f_bp(3);
    [~, dind] = min(abs(freq - bp_init));
end


%% Step 4. perform initial detrending
Y_init = Y(ind0);
f_init = freq(ind0);

% definie initial parameters and constrains
x0 = [logp(1), exp(logp(2)), logp(1), exp(logp(2)), bp_init];
lb = [-10,-inf,-10,-inf,bp_lo];
ub = [10,inf,10,inf,bp_hi];

% without weighted loss
% cost_init = @(x)sum((Y_init - (log(x(2)*f_init.^x(1)).*(f_init<=x(5)) + log(x(4)*f_init.^x(3)).*(f_init>x(5)))).^2);

% with weighted loss
% cost_init = @(x)sum((((Y_init - log(x(2)*f_init.^x(1))).*(f_init<=x(5))).^2)*(1-x(5)/range(f_init))) + sum((((Y_init - log(x(4)*f_init.^x(3))).*(f_init>x(5))).^2)*(x(5)/range(f_init)));
% cost_init = @(x)sum(((Y_init - log(x(2)*f_init.^x(1))).*(length(f_init)/sum(f_init<=x(5))).*(f_init<=x(5))).^2) +...
%     sum(((Y_init - log(x(4)*f_init.^x(3))).*(length(f_init)/sum(f_init>=x(5))).*(f_init>x(5))).^2);

% adjust weights for frequency distribution
W_freq_init = f_init.^(-m_set.f_w);
% W_freq_init = ones(size(f_init));

cost_init = @(x)sum(((Y_init - log(x(2)*f_init.^x(1))).*W_freq_init.*(f_init<=x(5))).^2) +...
    sum(((Y_init - log(x(4)*f_init.^x(3))).*W_freq_init.*(f_init>x(5))).^2);

b0 = fmincon(cost_init,x0,[],[],[],[],lb,ub,@bimodal_con,f_options);
A_init = log(b0(2).*freq.^b0(1)).*(freq<=b0(5))+log(b0(4).*freq.^b0(3)).*(freq>b0(5));

%% Step 5. detrend signal
Yd = Y - A_init;

% moving average smoothing
Yd = movmean(Yd,(1/f_res)*m_set.f_ma);

%% Step 6. find Gaussian peaks

% obtain initial parameters for Gaussian peaks
gfits = func_Gaussian_peaks(Yd, freq, m_set);

if ~isempty(gfits) % peaks present

    % drop redundant/extreme peaks
    gfits = func_drop_gfits(freq, gfits, m_set.th_edge, m_set.th_overlap, m_set.max_peak_N);

    if ~ isempty(gfits)

        % estimate multi-Gaussian fit
        [gfits_lb, gfits_ub] = func_gfits_constraints(gfits,m_set);

        gfits_multi = fmincon(@(x)func_multi_Gaussian_SSE(x,Yd,freq),gfits,[],[],[],[],gfits_lb,gfits_ub,[],f_options);
        gf_curve = func_multi_Gaussian(gfits_multi,freq);

        gpeaks = reshape(gfits_multi,[3,length(gfits_multi)/3])';
        gpeaks = sortrows(gpeaks,2,'ascend');

    else % no peaks
        gfits_multi = [];
        gf_curve = zeros(size(freq));
        gpeaks = [];
    end


else % no peaks
    gfits_multi = [];
    gf_curve = zeros(size(freq));
    gpeaks = [];
end

%% Step 7. remove Gaussian peaks to obtain broadband fractal spectrum
Yf = Y - gf_curve;

%% Step 8. adjust loss function weights to account for peaks and logarithmic scaling

% obtain regions outside peaks
df_ind = zeros(size(freq));
if ~isempty(gfits_multi) % peaks present
    df_ind = zeros(size(Yf));
    N = length(gfits_multi)/3;
    for n = 1:N
        g_mu = gfits_multi((n-1)*3+2);
        g_std = gfits_multi((n-1)*3+3);

        df_ind(freq>g_mu-m_set.th_exclude*g_std & freq<g_mu+m_set.th_exclude*g_std) = 1;
    end
end

% adjust weights for peaks
W_peak = zeros(size(Yf));

W_peak(df_ind==1) = sum(df_ind)/length(Yf);
W_peak(df_ind==0) = (length(Yf)-sum(df_ind))/length(Yf);
W_peak = W_peak.^2;

% adjust weights for frequency distribution
W_freq = freq.^(-m_set.f_w);
% W_freq = ones(size(freq));


%% Step 8. perform bimodal fit
freq_np = freq;
Y_np = Yf;

% definie initial parameters and constrains
x0 = [logp(1), exp(logp(2)), logp(1), exp(logp(2)), bp_init];
lb = [-10,-inf,-10,-inf,bp_lo];
ub = [10,inf,10,inf,bp_hi];


% account for under/oversampling in low/high frequency regime

% without weighted loss
% cost_fit = @(x)sum((((Y_np - log(x(2)*freq_np.^x(1))).*(freq_np<=x(5))).^2)) + sum((((Y_np - log(x(4)*freq_np.^x(3))).*(freq_np>x(5))).^2));

% with weighted loss
% cost_fit = @(x)sum((((Y_np - log(x(2)*freq_np.^x(1))).*(freq_np<=x(5))).^2)*(1-x(5)/range(freq))) + sum((((Y_np - log(x(4)*freq_np.^x(3))).*(freq_np>x(5))).^2)*(x(5)/range(freq)));
% cost_fit = @(x)sum(abs((Y_np - log(x(2)*freq_np.^x(1))).*(freq_np<=x(5))))*(1-x(5)/range(freq)) + sum(abs((Y_np - log(x(4)*freq_np.^x(3))).*(freq_np>x(5)))*(x(5)/range(freq)));
% cost_fit = @(x)sum(((Y_np - log(x(2)*freq_np.^x(1))).*(length(freq_np)/sum(freq_np<=x(5))).*(freq_np<=x(5))).^2) +...
%     sum(((Y_np - log(x(4)*freq_np.^x(3))).*(length(freq_np)/sum(freq_np>=x(5))).*(freq_np>x(5))).^2);
% cost_fit = @(x)mean(((Y_np - log(x(2)*freq_np.^x(1))).*(freq_np<=x(5))).^2) +...
%     mean(((Y_np - log(x(4)*freq_np.^x(3))).*(freq_np>x(5))).^2);
cost_fit = @(x)sum(((Y_np - log(x(2)*freq_np.^x(1))).*W_peak.*W_freq.*(freq_np<=x(5))).^2) +...
    sum(((Y_np - log(x(4)*freq_np.^x(3))).*W_peak.*W_freq.*(freq_np>x(5))).^2);

b = fmincon(cost_fit,x0,[],[],[],[],lb,ub,@bimodal_con,f_options);
A = log(b(2).*freq.^b(1)).*(freq<=b(5))+log(b(4).*freq.^b(3)).*(freq>b(5));

% compute unimodal fit
% b_uni = fit(freq_np,Y_np,'log');
% A_uni = b_uni(freq);
p_uni = exp(polyfit(log(freq_np),Y_np,1));
b_uni_tmp = fminunc(@(x)sum((Y_np - log(x(2)*freq_np.^(x(1)))).^2), [p_uni(1), p_uni(2)],optimoptions('fminunc','Display','off'));
b_uni = struct('a', b_uni_tmp(1), 'b', b_uni_tmp(2));
A_uni = log(b_uni.b.*freq.^(b_uni.a));

%% Step 9. confirm spectral bimodality

[h, p, stat] = func_bimodal_slope_test(freq_np, Y_np, b, m_set.bimodal_alpha);

if h == 1
    conc = 'bimodal';
else
    conc = 'unimodal';
end


%% Step X. reconstruct combined fractal + oscillatory fit
Y_mmspm = A + gf_curve;
Y_umspm = A_uni + gf_curve;

%% Step XX. store results



b_struct = struct(...
    'beta_lo', b(1),...
    'intercept_lo', b(2),...
    'beta_hi', b(3),...
    'intercept_hi', b(4),...
    'breakpoint', b(5),...
    'stat', stat,...
    'conclusion', conc,...
    'p', p,...
    'beta_uni', b_uni.a,...
    'intercept_uni', b_uni.b,...
    'freq', freq,...
    'Y_raw', Y_raw,...
    'Y', Y,...
    'A', A,...
    'A_uni', A_uni,...
    'peaks', gpeaks,...
    'Y_mmspm', Y_mmspm,...
    'Y_umspm', Y_umspm);

% store data for plotting
plot_struct = struct(...
    'freq', freq,...
    'Y_raw', Y_raw,...
    'Y', Y,...
    'E_lo', E_lo,...
    'Y_logfit', Y_logfit,...
    'dind', dind,...
    'bp_init', bp_init,...
    'ind0', ind0,...
    'Yd0', Yd0,...
    'f_init', f_init,...
    'Y_init', Y_init,...
    'A_init', A_init,...
    'b0', b0,...
    'Yd', Yd,...
    'gf_curve', gf_curve,...
    'Yf', Yf,...
    'Y_mmspm', Y_mmspm,...
    'A', A,...
    'b', b);

b_struct.plot_struct = plot_struct;


%% plot (debug)
if m_set.flag_fig
    plot_mmspm(b_struct);
end

%% display results
if m_set.flag_disp
    disp_mmspm(b_struct)
end

% disp(b0)
end

%% helper functions
function [c, ceq] = bimodal_con(x)
c = [];
ceq = x(2)*x(5)^x(1) - (x(4)*x(5)^x(3));
end
