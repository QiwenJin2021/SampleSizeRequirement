function [PSD_simu,D_simu] = calc_PSD_D_simu(DD,dm,dn,xmin,xmax)
%
% pdf = [fn,fa,fm];
% cdf = [Fn,Fa,Fm];

Pnum = length(DD);
binNum = length(dm);
pdf = zeros(binNum,3);
cdf = zeros(binNum,3);

DD1 = DD;
DD2 = DD.^2;
DD3 = DD.^3;
nSum = sum(DD1);
aSum = sum(DD2);
vSum = sum(DD3);
Dn = nSum/Pnum;
%Da = (aSum/Pnum)^0.5;
Dm = (vSum/Pnum)^(1/3);
D32 = vSum/aSum;

tmp = zeros(binNum,3);

ldm = [xmin; dm(1:end-1)];
udm = dm;

ldn = [xmin; dn(1:end-1)];
udn = dn;

for i = 1:Pnum
    for j = 1:binNum
    if ldm(j)<=DD(i) && udm(j)>DD(i)
        tmp(j,2) = tmp(j,2) + DD2(i);
        tmp(j,3) = tmp(j,3) + DD3(i);
    end
    if ldn(j)<=DD(i) && udn(j)>DD(i)
        tmp(j,1) = tmp(j,1) + 1;
    end

    end
end
pdf(:,1) = tmp(:,1)./Pnum;
pdf(:,2) = tmp(:,2)./aSum;
pdf(:,3) = tmp(:,3)./vSum;

cdf(1,:) = pdf(1,:);
for j = 2:binNum
    cdf(j,:) = cdf(j-1,:) + pdf(j,:);
end

fm = pdf(:,3);
fn = pdf(:,1);
Fm = cdf(:,3);
Fn = cdf(:,1);

PSD_simu = [Fm fm Fn fn];

[C, ia, ic] = unique(Fm);
x = dm(ia);
x=x';
F = griddedInterpolant(C, x, 'pchip');
dm10 = F(0.1);
dm50 = F(0.5);
dm90 = F(0.9);
Dm50 = [dm10 dm50 dm90];

[C, ia, ic] = unique(Fn);
x = dn(ia);
x=x';
F = griddedInterpolant(C, x, 'pchip');
dn10 = F(0.1);
dn50 = F(0.5);
dn90 = F(0.9);
Dn50 = [dn10 dn50 dn90];

% D_simu = [Dm50 Dn50 Dm Dn D32];
D_simu = [dm50 dn50 Dm Dn D32];

end

