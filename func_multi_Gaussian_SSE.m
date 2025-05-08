function [ SSE ] = func_multi_Gaussian_SSE(x, Y, freq)

% loss function (SSE) for multi-Gaussian fit
N = length(x)/3;

gcurve = zeros(size(Y));

for n = 1:N
    gf = x((n-1)*3+1:n*3);
    gcurve_n = gf(1)*exp(-(1/2)*(((freq-gf(2))/(gf(3))).^2));
    gcurve = gcurve + gcurve_n;
end

SSE = sum((Y - gcurve).^2);

end