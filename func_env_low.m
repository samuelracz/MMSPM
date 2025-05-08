function [ ylower ] = func_env_low(Y,W)

Y = Y(:);
W = W(:);

N = length(Y);
xt = (1:N)';
Nw = length(W);
ylower = zeros(N,Nw);

for w = 1:Nw
    wi = W(w);
    min_peaks = zeros(floor(N/wi),2);
    min_locs = zeros(floor(N/wi),2);

    % forward pass
    for i = 1:1:floor(N/wi)
        tstart = (i-1)*wi+1;
        tstop = i*wi;
        Xt = Y(tstart:tstop);
        [minval,minind] = min(Xt);
        min_peaks(i,1) = minval;
        min_locs(i,1) = minind + tstart-1;
    end

    % backward pass
    for i = 1:1:floor(N/wi)
        tstart = N-i*wi+1;
        tstop = N-(i-1)*wi;
        Xt = Y(tstart:tstop);
        [minval,minind] = min(Xt);
        min_peaks(floor(N/wi)-i+1,2) = minval;
        min_locs(floor(N/wi)-i+1,2) = minind + tstart-1;
    end

    % eliminate duplicates
    tab_min = sortrows(table([min_peaks(:,1); min_peaks(:,2)],[min_locs(:,1); min_locs(:,2)],'VariableNames',{'val','ind'}),'ind','ascend');
    [~, ia] = unique(tab_min.ind);
    tab_min = tab_min(ia,:);

    % add first element - assume local linearity
    if tab_min.ind(1) ~= 1
        y_seg = Y(1:W(1));
        x_seg = (0:W(1)-1)';
        yp = polyfit(x_seg,y_seg,1);

        % compute detrended variance
        y_std = std(y_seg - (yp(1)*x_seg + yp(2)));
        tab_min_tmp = [struct2table(struct('val', yp(2)-2*y_std ,'ind', 1)); tab_min];
        tab_min = tab_min_tmp;
    end

    % add last element - assume local linearity
    if tab_min.ind(end) ~= N
        y_seg = Y(end-W(1)+1:end);
        x_seg = (-W(1)+1:0)';
        yp = polyfit(x_seg,y_seg,1);

        % compute detrended variance
        y_std = std(y_seg - (yp(1)*x_seg + yp(2)));
        tab_min_tmp = [tab_min; struct2table(struct('val', yp(2)-2*y_std, 'ind', N))];
        tab_min = tab_min_tmp;
    end

    % interpolation
    ylower(:,w) = makima(tab_min.ind,tab_min.val,xt);
end

end