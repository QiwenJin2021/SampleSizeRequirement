function [X] = gen_sample(method,type,N,n,xmin,xmax,De)
% Generate particle samples X
% N: total particle number
% Methods: Inverse Transform Method (ITM) or Acceptance-Rejection Method (ARM)
% Type: GSD2 and RRD2 are modified GSD and RRD with xmin and xmax.
switch method
    case 'ITM'
        u = rand(1,N);
        u=u';
        switch type
            case 'GSD'
                X = xmax*(u).^(1/(n-3));
            case 'PD'
                X = xmin./(1-u).^(1/(n+3));
            case 'RRD'
                xn = gammaincinv(1-((1-u)*igamma(1-3/n,0)/gamma(1-3/n)),1-3/n,'upper');
                X = xn.^(1/n)*De;
            case 'RRD2'
                aa = 1-3/n;
                bb = (xmin/De)^n;
                cc = (xmax/De)^n;
                t1 = igamma(aa,bb)-(igamma(aa,bb)-igamma(aa,cc))*u;
                xn = gammaincinv(t1/gamma(aa),aa,'upper');
                X = xn.^(1/n)*De;
            case 'GSD2'
                H = hypergeom([3,n],1+n,(xmin-xmax)/xmin);

                for i = 1:N
                    syms x
                    eqn = u(i) == (x-xmin)^n * hypergeom([3,n],1+n,(xmin-x)/xmin)/(xmax-xmin)^n / H;
                    s = vpasolve(eqn, x, [xmin,xmax]);
                    X(i) = double(s);
                end
        end

    case 'ARM'
        switch type
            case 'GSD2'
                notenough = 1;
                A = xmin^3 / hypergeom([3,n],1+n,(xmin-xmax)/xmin);
%                 x = xmin : 0.1 : xmax;
%                 y = (x-xmin).^(n-1).*x.^(-3)*n/(xmax-xmin)^n * A;
%                 ymax = max(y);
                fun = @(x) -1 * A * n / (xmax-xmin) * ((x-xmin)/(xmax-xmin))^(n-1) * x^(-3);
                [x,ymax] = fminbnd(fun,xmin,xmax);
                ymax = -1 * ymax;
                % 舍选法
                X1 = [];
                while notenough
                    ux = rand(1,2*N);
                    XX = ux*(xmax-xmin)+xmin;
                    uy = rand(1,2*N);
                    uy = uy*(ymax+0.0001);
                    YY = (XX-xmin).^(n-1).*XX.^(-3)*n/(xmax-xmin)^n * A;
                    indX = find(uy <= YY);
                    X2 = XX(indX);
                    X1 = [X1 X2];

                    if length(X1)>=N
                        notenough = 0;
                    end

                end
                X = datasample(X1,N);

            case 'RRD2'
                notenough = 1;

                e1 = exp(-(xmin/De)^n);
                e2 = exp(-(xmax/De)^n);
                 g1 = igamma(1-3/n,(xmin/De)^n);
                 g2 = igamma(1-3/n,(xmax/De)^n );

                g1 = real(g1);
                g2 = real(g2);
                A = ( e1-e2)/( g1-g2 );

%                 x = xmin : 0.01 : xmax;
%                 e = exp(-(x/De).^n);
%                 y = n / De^n * x.^(n-4) .* e /(e1-e2) * A;
%                 ymax = 1.2*max(y);

                %fun = @(x) n / De^n * x^(n-4) * exp(-(x/De).^n) /(e1-e2) * A;
                fun = @(x) -1 * n / De^(n-3) * x^(n-4) * exp(-(x/De).^n) /(g1-g2);
                [x,ymax] = fminbnd(fun,xmin,xmax);
                ymax = -1 * ymax;

                % 舍选法
                X1 = [];
                while notenough
                    ux = rand(1,2*N);
                    XX = ux*(xmax-xmin)+xmin;
                    uy = rand(1,2*N);
                    uy = uy*(ymax+0.0001);
                    E = exp(-(XX/De).^n);
                    %YY = n / De^n * XX.^(n-4) .* E /(e1-e2) * A;
                    YY = n / De^(n-3) * XX.^(n-4) .* E /(e1-e2) * A;
                    indX = find(uy <= YY);
                    X2 = XX(indX);
                    X1 = [X1 X2];

                    if length(X1)>=N
                        notenough = 0;
                    end

                end
                X = datasample(X1,N);


        end
end

X = X';

end
