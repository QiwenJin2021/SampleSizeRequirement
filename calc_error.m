function [err_P,err_avg] = calc_error(errs,P)
%UNTITLED6 此处提供此函数的摘要

[N1, N2]= size(errs);
errs1 = abs(errs);

err_max = max(max(errs1));
% delta_N = 200;
% delta = err_max/delta_N;
delta = 0.0001;
err_max = ceil(err_max/delta)*delta;
bins = err_max / delta;
bins = ceil(bins);
% bins = delta_N;
err_delta = delta:delta:err_max;

err_avg = sum(errs1,1)/N1;

err_dis = zeros(bins,N2);
err_dis_norm = zeros(bins,N2);
err_dis_cum = zeros(bins,N2);


err_P = zeros(1,N2);

for i = 1:N2
    for j = 1:N1
        K = errs1(j,i) / delta;
        K = ceil(K);
        if K==0
            err_dis(1,i) = err_dis(1,i) + 1;
        end
        if K==bins
            err_dis(bins,i) = err_dis(bins,i) + 1;
        end
        if K>0 && K<bins
            err_dis(K,i) = err_dis(K,i) + 1;
        end
    end
    err_dis_norm(:,i) = err_dis(:,i) / N1;
    err_dis_cum(1,i) = err_dis_norm(1,i);
    for kk = 2:bins
        err_dis_cum(kk,i) = err_dis_cum(kk-1,i) + err_dis_norm(kk,i);
    end
    [C, ia, ic] = unique(err_dis_cum(:,i));
    x = err_delta(ia);
    x=x';
    F = griddedInterpolant(C, x, 'pchip');
    err_P(1,i) = F(P);

end