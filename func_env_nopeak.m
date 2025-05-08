function [inds_out] = func_env_nopeak(Y,freq,W,th_height,th_width)

pmod = polyfit(log(freq),Y,1);
Yd = pmod(1)*log(freq)+pmod(2);

[pks, locs, ws, ~] = findpeaks(movmean(Y-Yd,10*W),freq,'MinPeakProminence',th_height,'MinPeakWidth',W,'Annotate','extents','WidthReference','halfheight');

if isempty(pks) % no peaks
    inds_out = [];
else
    inds = ones(length(Y),1);

    for p = 1:length(pks)
        inds(freq >= locs(p)-th_width*ws(p) & freq <= locs(p)+th_width*ws(p)) = 0;
    end

    inds_out = inds == 1;

    % for safety, always include last elements
    inds_out(1:4) = 1;
    inds_out(end-4:end) = 1;
end


% plot(freq,Y)
% hold on
% plot(freq(inds==1),Y(inds==1),'r*')
% plot(locs,pks,'go')

end