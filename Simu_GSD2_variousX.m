% Simulation of GSD
% 回复审稿人意见补充数据
% 1、验证生成大样本的准确性
% 2、不同模拟次数的影响
% 3、不同定义方式产生的跨度值区别，数目分布的跨度值先增后减的原因
%%%% 4、不同粒径边界值的影响，是否有普适性  %%%%
% 5、Gy采样理论计算结果

% 针对第四点：不同粒径边界值的影响，是否有普适性
clear all;clc;
addpath(genpath(pwd));


% 4. 不同xmin、xmax、De下，普适性
type = 'GSD2';
method = 'ARM';
n = 3.2;
xmin = 1;
xmax = 20:20:160;
De = 50;

P = 0.95;
filepath = ['D:\Mycode\ParticleAnalysis\PSD_simu_v2\SampleData\', type, '_', method, '_bounds'  '\'];

for kk = 1:length(xmax)
    SR1(kk,:) = calc_SR(type,n,xmin,xmax(kk),De);
%     SR2(kk,:) = calc_SR2(type,n,xmin,xmax(kk),De);
%     SR3(kk,:) = calc_SR3(type,n,xmin,xmax(kk),De);

end

N = 10000000;

for j = 1:length(xmax)
    fprintf('xmax=%2d\n', xmax(j));
    X = gen_sample(method,type,N,n,xmin,xmax(j),De);
    fname = [filepath type '_n_'  num2str(n) '_' num2str(xmin) '_' num2str(xmax(j))  '.mat'];
    save(fname, 'X' );

end


% xm1 = xmin:xmax-1;
% xm1 = xm1';
% xn1 =xm1 +0.5 ;
% for kk = 1:length(xmax)
%     [xm,xn] = get_bin(type,n,xmin,xmax(kk),De);
%     PSD_the = calc_PSD_theory(type,xm,xn,n,xmin,xmax(kk),De);
%     D_the = calc_D_theory(type,n,xmin,xmax(kk),De);
%     fname = [filepath type '_n_'  num2str(n) '_' num2str(xmin) '_' num2str(xmax(kk))  '.mat'];
%     load(fname);
%     [PSD_simu,D_simu(kk,:)] = calc_PSD_D_simu(X,xm,xn,xmin,xmax(kk));
%     Fm_the(:,kk) = PSD_the(:,1);
%     Fn_the(:,kk) = PSD_the(:,3);
%     Fm_sim(:,kk) = PSD_simu(:,1);
%     Fn_sim(:,kk) = PSD_simu(:,3);
% 
%     fm_the(:,kk) = PSD_the(:,2);
%     fn_the(:,kk) = PSD_the(:,4);
% 
%     fm_sim(:,kk) = PSD_simu(:,2);
%     fn_sim(:,kk) = PSD_simu(:,4);
% end



% SampleSize = [1000:1000:20000];
SampleSize = 10000;
SampleSize = SampleSize';
len_SampleSize = length(SampleSize);

times = 2000;
for kk = 1:length(xmax)
    [xm,xn] = get_bin(type,n,xmin,xmax(kk),De);
    D_the = calc_D_theory(type,n,xmin,xmax(kk),De);
    fname = [filepath type '_n_'  num2str(n) '_' num2str(xmin) '_' num2str(xmax(kk)) '.mat'];
    load(fname);
    fprintf('xmax=%2d\n', xmax(kk));
%     PSD_the = calc_PSD_theory(type,xm1,xn1,n,xmin,xmax(kk),De(kk));
    for j = 1:times
        %             fprintf('n=%2d,i=%2d,j=%2d\n', n(kk),i,j);
        X1 = datasample(X,SampleSize);
        [PSD_simu,D_simu(j,:)] = calc_PSD_D_simu(X1,xm,xn,xmin,xmax(kk));
        err_D(j,:) = D_simu(j,:) ./ D_the - 1;
    end
    %         err_D_abs = abs(err_D);
    [err_P(kk,:),err_avg(kk,:)] = calc_error(err_D,P);
    err_Dm50(:,kk) = err_D(:,1);
    err_Dn50(:,kk) = err_D(:,2);
    err_Dm(:,kk) = err_D(:,3);
    err_Dn(:,kk) = err_D(:,4);
    err_D32(:,kk) = err_D(:,5);

    %     matfile = [filepath 'Results_GSD2_n_' num2str(n) 'xmax_' num2str(xmax(kk)) '.mat'];
    %     save(matfile,'err_Dm50','err_Dn50','err_Dm','err_Dn','err_D32','err_P','err_avg');
end