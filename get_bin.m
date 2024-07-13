function [dm,dn] = get_bin(type,n,xmin,xmax,De)
% Calculate the x-coordinate array of PSD.
% d: particle size.
% Fm: mass-based cumulative distribution function.
% GSD2 and RRD2 are modified GSD and RRD with xmin and xmax.

Fm = 0.01:0.01:0.99;
Fn = 0.01:0.01:0.99;
k = length(Fn);

switch type
    case 'GSD2'
        dm = Fm.^(1/n)*(xmax-xmin)+xmin;
        dm = [dm xmax];

        H = hypergeom([3,n],1+n,(xmin-xmax)/xmin);
        for i=1:k
            syms d1
            eqn = Fn(i) == (d1-xmin)^n * hypergeom([3,n],1+n,(xmin-d1)/xmin)/(xmax-xmin)^n / H;
            s = vpasolve(eqn,d1,[xmin,xmax]);
            dn(i) = double(s);
        end
        dn = [dn xmax];

    case 'PD'
        dm = xmin./(1-Fm).^(1/n);
        dmax1 = xmin /(1-0.9999).^(1/n);
        dm = [dm dmax1];

        dn = ((1-Fn).*xmin^(-(n+3))).^(-1/(n+3));
        dmax2 = ((1-0.9999)*xmin^(-(n+3)))^(-1/(n+3));
        dn = [dn dmax2];

    case 'RRD2'
        bb = exp(-(xmin/De)^n);
        cc = exp(-(xmax/De)^n);
        dm = De.*(-log(bb-Fm*(bb-cc))).^(1./n);
        dm = [dm xmax];

        aa = 1-3/n;
        bb = (xmin/De)^n;
        cc = (xmax/De)^n;
%         t1 = igamma(aa,bb)-(igamma(aa,bb)-igamma(aa,cc))*Fn;
%         xn = gammaincinv(t1/gamma(aa),aa,'upper');
%         dn = xn.^(1/n)*De;
%         dn = [dn xmax];

        for i=1:k
            syms d1
            eqn = Fn(i) == (igamma(aa,bb)-igamma(aa,(d1/De)^n))/(igamma(aa,bb)-igamma(aa,cc));
            s = vpasolve(eqn,d1,[xmin,xmax]);
            dn(i) = double(s);
        end

        dn = [dn xmax];

    case 'GSD'
        dm = Fm.^(1/n)*xmax;
        dm = [dm xmax];

        dn = Fn.^(1/(n-3)) * xmax;
        dn = [dn xmax];

    case 'RRD'
        dm = De.*(-log(1-Fm)).^(1./n);
        dmax1 = De.*(-log(1-0.9999)).^(1./n);
        dm = [dm dmax1];

        dn = gammaincinv(1-((1-Fn)*igamma(1-3/n,0)/gamma(1-3/n)),1-3/n,'upper');
        dmax2 = gammaincinv(1-((1-0.9999)*igamma(1-3/n,0)/gamma(1-3/n)),1-3/n,'upper');
        dn = [dn dmax2];
end
dm =dm';
dn = dn';

end