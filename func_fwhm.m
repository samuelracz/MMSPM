function [fm, fwhm_neg, fwhm_pos] = func_fwhm(xt,Y)

[fm, max_ind] = max(Y);

X_neg = Y(1:max_ind);
X_pos = Y(max_ind:end);

dX_neg = abs(X_neg-fm/2);
[~,neg_ind] = min(dX_neg);

dX_pos = abs(X_pos-fm/2);
[~,pos_ind] = min(dX_pos);

fwhm_neg = xt(max_ind) - xt(neg_ind);
fwhm_pos = xt(max_ind + pos_ind) - xt(max_ind);

end