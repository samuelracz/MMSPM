function [] = func_env_skew(Y,freq,W)

Y = log(Y(:));
freq = log(freq(:));
W = W(:);

N = length(Y);
Nw = length(W);

yskew = zeros(N,Nw);
ylower = zeros(N,Nw);

for w = 1:Nw
    ns = floor(N/W(w));
    tmp_fp = zeros(N,1);
    tmp_bp = zeros(N,1);

    % forward pass
    for s = 1:ns
        Y_tmp = Y((s-1)*W(w)+1:s*W(w));
        Y_dt = Y_tmp - linspace(Y_tmp(1),Y_tmp(end),W(w))';
        tmp_fp((s-1)*W(w)+1:s*W(w)) = abs(sum(Y_dt));
    end

    if N/W(w) > ns*W(w)
        Y_tmp = Y(ns*W(w)+1:end);
        Y_dt = Y_tmp - linspace(Y_tmp(1),Y_tmp(end),W(w))';
        tmp_fp(ns*W(w)+1:end) = abs(sum(Y_dt));
    end

    % backward pass
    Y_flip = flipud(Y);
    for s = ns
        Y_tmp = Y_flip((s-1)*W(w)+1:s*W(w));
        Y_dt = Y_tmp - linspace(Y_tmp(1),Y_tmp(end),W(w))';
        tmp_bp((s-1)*W(w)+1:s*W(w)) = abs(sum(Y_dt));
    end

    if N/W(w) > ns*W(w)
        Y_tmp = Y_flip(ns*W(w)+1:end);
        Y_dt = Y_tmp - linspace(Y_tmp(1),Y_tmp(end),W(w))';
        tmp_bp(ns*W(w)+1:end) = abs(sum(Y_dt));
    end
    tmp_bp = flipud(tmp_bp);

    yskew(:,w) = mean([tmp_fp, tmp_bp],2);
end

subplot(211)
plot(freq,Y)
subplot(212)
plot(freq,yskew)

end