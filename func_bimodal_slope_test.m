function [h, p, stat] = func_bimodal_slope_test(freq_np, Y_np, b, test_alpha)

%% obtain main variables

% get fit data for the two ranges
f_lo = freq_np(freq_np <= b(5));
Y_lo = Y_np(freq_np <= b(5));
A_lo = log(b(2)*f_lo.^(b(1)));
N_lo = length(f_lo);

f_hi = freq_np(freq_np > b(5));
Y_hi = Y_np(freq_np > b(5));
A_hi = log(b(4)*f_hi.^(b(3)));
N_hi = length(f_hi);

% get residuals
e_lo = Y_lo - A_lo;
e_hi = Y_hi - A_hi;

% get residual sum of squares
RSS_lo = sum(e_lo.^2);
RSS_hi = sum(e_hi.^2);

% get variance of errors
sigma_lo = RSS_lo/(N_lo-2);
sigma_hi = RSS_hi/(N_hi-2);

% get variance of data
var_lo = sum((Y_lo - mean(Y_lo)).^2)/N_lo;
var_hi = sum((Y_hi - mean(Y_hi)).^2)/N_hi;

% get standard error of the slope
SE_lo = sqrt(sigma_lo/var_lo);
SE_hi = sqrt(sigma_hi/var_hi);

%% hypothesis testing

% standard error of slope difference
SE_bimodal = sqrt(SE_lo.^2 + SE_hi.^2);

% generating test statistic
t_bimodal = (b(1) - b(3))/SE_bimodal;

% degrees of freedom
df_bimodal = min([N_lo-2,N_hi-2]);

% cumulative distribution function
tcdf_bimodal = tcdf(t_bimodal,df_bimodal);

% get p-value
if t_bimodal >= 0
    p_bimodal = 1-tcdf_bimodal;
elseif t_bimodal < 0
    p_bimodal = tcdf_bimodal;
end

% make decision
if p_bimodal < test_alpha
    h_bimodal = 1;
else
    h_bimodal = 0;
end

%% store output

h = h_bimodal;
p = p_bimodal;
stat = struct(...
    'p', p,...
    'h', h,...
    'tstat', t_bimodal,...
    'df', df_bimodal,...
    'SE', SE_bimodal);

% figure;
% plot(log(freq_np),Y_np)
% hold on
% plot(log(f_lo),Y_lo,'*')
% plot(log(f_hi),Y_hi,'.')
% plot(log(f_lo),A_lo,'r-','LineWidth',2)
% plot(log(f_hi),A_hi,'g-','LineWidth',2)


end