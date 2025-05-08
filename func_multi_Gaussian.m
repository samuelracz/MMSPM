function [ gcurve ] = func_multi_Gaussian(gfits, freq)

% generate multi-Gaussian curve
N = length(gfits)/3;

gcurve = zeros(size(freq));

for n = 1:N
    gf = gfits((n-1)*3+1:3*n);
    gcurve_n = gf(1)*exp(-(1/2)*(((freq-gf(2))/(gf(3))).^2));
    gcurve = gcurve + gcurve_n;
end

end