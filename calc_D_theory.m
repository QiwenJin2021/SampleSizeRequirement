function [D_the] = calc_D_theory(type,n,xmin,xmax,De)
%

u = [0.1 0.5 0.9];
switch type
    % GSD with xmin and xmax
    case 'GSD2'
        % mass-weighted percentile diameter
        Dm50 = u.^(1/n)*(xmax-xmin)+xmin;

        H = hypergeom([3,n],1+n,(xmin-xmax)/xmin);
        A = xmin^3 / H;

        syms d1 d2 d3
        eqn1 = 0.1 == (d1-xmin)^n * hypergeom([3,n],1+n,(xmin-d1)/xmin)/(xmax-xmin)^n / H;
        eqn2 = 0.5 == (d2-xmin)^n * hypergeom([3,n],1+n,(xmin-d2)/xmin)/(xmax-xmin)^n / H;
        eqn3 = 0.9 == (d3-xmin)^n * hypergeom([3,n],1+n,(xmin-d3)/xmin)/(xmax-xmin)^n / H;
        dn10 = vpasolve(eqn1,d1,[xmin,xmax]);
        dn50 = vpasolve(eqn2,d2,[xmin,xmax]);
        dn90 = vpasolve(eqn3,d3,[xmin,xmax]);
        dn10 = double(dn10);
        dn50 = double(dn50);
        dn90 = double(dn90);

        Dn50 = [dn10 dn50 dn90];

        fun3 = @(x) x.^3 * A * n / (xmax-xmin) .* ((x-xmin)/(xmax-xmin)).^(n-1) .* x.^(-3);
        q3 = integral(fun3,xmin,xmax);
        Dm = q3^(1/3);

        fun1 = @(x) x .* A * n / (xmax-xmin) .* ((x-xmin)/(xmax-xmin)).^(n-1) .* x.^(-3);
        q1 = integral(fun1,xmin,xmax);
        Dn = q1;

        fun2 = @(x) x.^2 * A * n / (xmax-xmin) .* ((x-xmin)/(xmax-xmin)).^(n-1) .* x.^(-3);
        q2 = integral(fun2,xmin,xmax);
        D32 = q3 / q2;



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
        eqn1 = 0.1 == (g1-igamma(aa,(d1/De)^n))/(g1-g2);
        eqn2 = 0.5 == (g1-igamma(aa,(d2/De)^n))/(g1-g2);
        eqn3 = 0.9 == (g1-igamma(aa,(d3/De)^n))/(g1-g2);
        dn10 = vpasolve(eqn1,d1,[xmin,xmax]);
        dn50 = vpasolve(eqn2,d2,[xmin,xmax]);
        dn90 = vpasolve(eqn3,d3,[xmin,xmax]);
        dn10 = double(dn10);
        dn50 = double(dn50);
        dn90 = double(dn90);
        Dn50 = [dn10 dn50 dn90];

        fun3 = @(x) x.^3 * A * n / De^n .* x.^(n-4) .* exp(-(x/De).^n) / (e1 -e2);
        q3 = integral(fun3,xmin,xmax);
        Dm = q3^(1/3);

        fun1 = @(x) x .* A * n / De^n .* x.^(n-4) .* exp(-(x/De).^n) / (e1 -e2);
        q1 = integral(fun1,xmin,xmax);
        Dn = q1;

        fun2 = @(x) x.^2 * A * n / De^n .* x.^(n-4) .* exp(-(x/De).^n) / (e1 -e2);
        q2 = integral(fun2,xmin,xmax);
        D32 = q3 / q2;

        % PD
    case 'PD'

        Dm50 = xmin./(1-u).^(1/n);

        Dn50 = ((1-u).*xmin^(-(n+3))).^(-1/(n+3));

        A = (n+3)*xmin^3/n;

        fun3 = @(x) x.^3 .* A * n * xmin^n ./ x.^(n + 4);
        q3 = integral(fun3,xmin,Inf);
        Dm = q3^(1/3);

        fun1 = @(x) x .* A * n * xmin^n ./ x.^(n + 4);
        q1 = integral(fun1,xmin,Inf);
        Dn = q1;

        fun2 = @(x) x.^2 .* A * n * xmin^n ./ x.^(n + 4);
        q2 = integral(fun2,xmin,Inf);
        D32 = q3 / q2;

        % GSD
    case 'GSD'

        Dm50 = u.^(1/n)*xmax;

        Dn50 = u.^(1/(n-3)) * xmax;

        A = xmax^3 * (n - 3) / n;

        fun3 = @(x) x.^3 .* A * n / xmax^n .* x^(n - 4);
        q3 = integral(fun3,0,xmax);
        Dm = q3^(1/3);

        fun1 = @(x) x .* A * n / xmax^n .* x^(n - 4);
        q1 = integral(fun1,0,xmax);
        Dn = q1;

        fun2 = @(x) x.^2 .* A * n / xmax^n .* x^(n - 4);
        q2 = integral(fun2,0,xmax);
        D32 = q3 / q2;

        % RRD
    case 'RRD'

        Dm50 = De*(-log(1-u)).^(1/n);

        Dn50 = gammaincinv(1-((1-u)*igamma(1-3/n,0)/gamma(1-3/n)),1-3/n,'upper');

        A = De^3 / igamma(1-3/n,0);

        fun3 = @(x) x.^3 .* A * n / De^n .* x^(n - 4) .* exp(-(x/De).^n);
        q3 = integral(fun3,0,Inf);
        Dm = q3^(1/3);

        fun1 = @(x) x .* A * n / De^n .* x^(n - 4) .* exp(-(x/De).^n);
        q1 = integral(fun1,0,Inf);
        Dn = q1;

        fun2 = @(x) x.^2 .* A * n / De^n .* x^(n - 4) .* exp(-(x/De).^n);
        q2 = integral(fun2,0,Inf);
        D32 = q3 / q2;


end

% D_the = [Dm50 Dn50 Dm Dn D32];
D_the = [Dm50(2) Dn50(2) Dm Dn D32];

end