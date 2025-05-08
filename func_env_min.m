function [ ylower ] = func_env_min(Y, freq, m_set)

N = length(Y);
fmins = [];
ymins = [];

% find rolling minimum
for n = 1:N
    if n == 1
        fmins = freq(n);
        ymins = Y(n);
    else
        if Y(n) < ymins(end)
            ymins_tmp = [ymins; Y(n)];
            ymins = ymins_tmp;

            fmins_tmp = [fmins; freq(n)];
            fmins = fmins_tmp;
        end
    end
end

% check edges
if fmins(end) ~= freq(end)
    ymins_tmp = [ymins; ymins(end)];
    ymins = ymins_tmp;

    fmins_tmp = [fmins; freq(end)];
    fmins = fmins_tmp;
end

% smooth curve
ylower = interp1(fmins,ymins,freq);
% ylower = movmean(ylower, ceil(m_set.f_ma/m_set.f_res));



% ylower = makima(fmins,ymins,freq);
% ylower = interp1(fmins,ymins,freq);
% ylower = pchip(fmins,ymins,freq);

% figure('color','w','units','normalized','outerposition',[0.1 0.3 0.8 0.5]);
% subplot(131)
% plot(log(freq),log(Y))
% hold on
% plot(log(freq),log(ylower))
% 
% subplot(132)
% plot(log(freq),log(Y))
% hold on
% plot(log(freq),log(ylower2))
% 
% subplot(133)
% plot(log(freq),log(Y))
% hold on
% plot(log(freq),log(ylower3))

end