function [ gfits_out ] = func_drop_gfits(freq,gfits,th_edge,th_overlap,Nmax)

gfits_out = gfits;
drop_ind = zeros(length(gfits),1);
N = length(gfits_out)/3;

% drop peaks too close to the edge
for n = 1:N
    g_mu = gfits_out((n-1)*3+2);
    g_std = gfits_out((n-1)*3+3);

    if (g_mu-th_edge*g_std < min(freq)) || (g_mu+th_edge*g_std > max(freq))
        drop_ind((n-1)*3+1:3*n) = 1;
    end
end

gfits_out(drop_ind==1) = [];

% if no peaks remain, return
if isempty(gfits_out)
    return
end

drop_ind = zeros(length(gfits),1);
N = length(gfits_out)/3;

% if only one peak remains, return
if N == 1
    return
end

% drop peaks that overlap too much
for n1 = 1:N
    g_amp1 = gfits_out((n1-1)*3+1);
    g_mu1 = gfits_out((n1-1)*3+2);
    g_std1 = gfits_out((n1-1)*3+3);
    for n2 = n1+1:N
        g_amp2 = gfits_out((n2-1)*3+1);
        g_mu2 = gfits_out((n2-1)*3+2);
        g_std2 = gfits_out((n2-1)*3+3);
        if (abs(g_mu1-g_mu2) < th_overlap*g_std1) || (abs(g_mu1-g_mu2) < th_overlap*g_std2)
            if g_amp1 < g_amp2
                drop_ind((n1-1)*3+1:(n1-1)*3+3) = 1;
            else
                drop_ind((n2-1)*3+1:(n2-1)*3+3) = 1;
            end
        end
    end
end

gfits_out(drop_ind==1) = [];
N = length(gfits_out)/3;

% keep only the N largest peaks
if N > Nmax
    amps = gfits_out(1:3:end);
    amps_include = maxk(amps,Nmax);
    drop_ind = zeros(length(gfits_out),1);
    for n = 1:N
        if ~ismember(gfits((n-1)*3+1),amps_include)
            drop_ind((n1-1)*3+1:(n1-1)*3+3) = 1;
        end
    end
    gfits_out(drop_ind==1) = [];
end


end