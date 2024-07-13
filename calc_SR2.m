function [SR] = calc_SR2(type,n,xmin,xmax,De)

% Define: SR = (D90-D10)/D50

u = [0.1 0.5 0.9];
switch type
    % GSD with xmin and xmax
    case 'GSD2'
        % mass-weighted percentile diameter
        Dm50 = u.^(1/n)*(xmax-xmin)+xmin;

        H = hypergeom([3,n],1+n,(xmin-xmax)/xmin);

        syms d1 d2 d3
        eqn1 = u(1) == (d1-xmin)^n * hypergeom([3,n],1+n,(xmin-d1)/xmin)/(xmax-xmin)^n / H;
        eqn2 = u(2) == (d2-xmin)^n * hypergeom([3,n],1+n,(xmin-d2)/xmin)/(xmax-xmin)^n / H;
        eqn3 = u(3) == (d3-xmin)^n * hypergeom([3,n],1+n,(xmin-d3)/xmin)/(xmax-xmin)^n / H;
        dn1 = vpasolve(eqn1,d1,[xmin,xmax]);
        dn2 = vpasolve(eqn2,d2,[xmin,xmax]);
        dn3 = vpasolve(eqn3,d3,[xmin,xmax]);
        dn1 = double(dn1);
        dn2 = double(dn2);
        dn3 = double(dn3);

        SRm = (Dm50(3)-Dm50(1))/Dm50(2);
        SRn = (dn3-dn1)/dn2;
        SR = [SRm SRn];

        % RRD with xmin and xmax
    case 'RRD2'

        e1 = exp(-(xmin/De)^n);
        e2 = exp(-(xmax/De)^n);

        Dm50 = De*(-1*log(e1-u.*(e1-e2))).^(1/n);

        %e = exp(-(x/De).^n);

%         g = igamma(1-3/n,(x/De).^n );
        g1 = igamma(1-3/n,(xmin/De)^n);
        g2 = igamma(1-3/n,(xmax/De)^n );
%         g = real(g);
        g1 = real(g1);
        g2 = real(g2);
        A = De^3 * ( e1-e2)/( g1-g2 );

        aa = 1-3/n;
        bb = (xmin/De)^n;
        cc = (xmax/De)^n;
%         u = [0.1 0.5 0.9];
%         t1 = igamma(aa,bb)-(igamma(aa,bb)-igamma(aa,cc))*u;
%         xn = gammaincinv(t1/gamma(aa),aa,'upper');
%         Dn50 = xn.^(1/n)*De;

        syms d1 d2 d3
        eqn1 = u(1) == (g1-igamma(aa,(d1/De)^n))/(g1-g2);
        eqn2 = u(2) == (g1-igamma(aa,(d2/De)^n))/(g1-g2);
        eqn3 = u(3) == (g1-igamma(aa,(d3/De)^n))/(g1-g2);
        dn1 = vpasolve(eqn1,d1,[xmin,xmax]);
        dn2 = vpasolve(eqn2,d2,[xmin,xmax]);
        dn3 = vpasolve(eqn3,d3,[xmin,xmax]);
        dn1 = double(dn1);
        dn2 = double(dn2);
        dn3 = double(dn3);

        SRm = (Dm50(3)-Dm50(1))/Dm50(2);
        SRn = (dn3-dn1)/dn2;
        SR = [SRm SRn];

       
        % PD
    case 'PD'

        Dm50 = xmin./(1-u).^(1/n);

        Dn50 = ((1-u).*xmin^(-(n+3))).^(-1/(n+3));

        SRm = (Dm50(3)-Dm50(1))/Dm50(2);
        SRn = (Dn50(3)-Dn50(1))/Dn50(2);
        SR = [SRm SRn];

        % GSD
    case 'GSD'

        Dm50 = u.^(1/n)*xmax;

        Dn50 = u.^(1/(n-3)) * xmax;

        SRm = (Dm50(3)-Dm50(1))/Dm50(2);
        SRn = (Dn50(3)-Dn50(1))/Dn50(2);
        SR = [SRm SRn];
        

        % RRD
    case 'RRD'

        Dm50 = De*(-log(1-u)).^(1/n);

        Dn50 = gammaincinv(1-((1-u)*igamma(1-3/n,0)/gamma(1-3/n)),1-3/n,'upper');

        SRm = (Dm50(3)-Dm50(1))/Dm50(2);
        SRn = (Dn50(3)-Dn50(1))/Dn50(2);
        SR = [SRm SRn];


end



end