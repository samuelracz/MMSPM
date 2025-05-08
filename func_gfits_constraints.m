function [ gfits_lb, gfits_ub ] = func_gfits_constraints(gfits, m_set)

N = length(gfits)/3;

gfits_lb = zeros(size(gfits));
gfits_ub = zeros(size(gfits));

for n = 1:N
    % lower boundaries
    gfits_lb((n-1)*3+1) = m_set.min_peak_Amp; % lower boundary on peak amplitude
    gfits_lb((n-1)*3+2) = gfits((n-1)*3+2)-m_set.th_dist*(gfits((n-1)*3+3)); % lower boundary on center of peak
    gfits_lb((n-1)*3+3) = 0; % no lower boundary on SD

    % upper boundaries
    gfits_ub((n-1)*3+1) = inf; % no upper boundary on peak amplitude
    gfits_ub((n-1)*3+2) = gfits((n-1)*3+2)+m_set.th_dist*(gfits((n-1)*3+3)); % upper boundary on center of peak
    gfits_ub((n-1)*3+3) = 2*gfits((n-1)*3+3); % upper boundary on SD
end

end