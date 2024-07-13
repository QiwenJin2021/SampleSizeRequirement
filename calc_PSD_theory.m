function [PSD_the] = calc_PSD_theory(type,xm,xn,n,xmin,xmax,De)
%

Fm = zeros(length(xm),1);
fm = zeros(length(xm),1);
Fn = zeros(length(xn),1);
fn = zeros(length(xn),1);

switch type
    % GSD with xmin and xmax
    case 'GSD2'
        A = xmin^3 / hypergeom([3,n],1+n,(xmin-xmax)/xmin);
        for i = 1:length(xm)
            if xm(i)>xmax
                Fm(i) = 1;
            end
            if xm(i)<=xmax
                Fm(i) = ((xm(i)-xmin)/(xmax-xmin)).^n;
            end

             if xn(i)>xmax
                Fn(i) = 1;
            end
            if xn(i)<=xmax             
                Fn(i) = (xn(i)-xmin).^n.*hypergeom([3,n],1+n,(xmin-xn(i))/xmin)/(xmax - xmin)^n / xmin^3 * A;
            end

        end

        fm(1) = Fm(1)-0;
        fn(1) = Fn(1)-0;
        for i = 2:length(xm)
            fm(i) = Fm(i) - Fm(i-1);
            fn(i) = Fn(i) - Fn(i-1);
        end

        % RRD with xmin and xmax
    case 'RRD2'
        em = exp(-(xm/De).^n);
        en = exp(-(xn/De).^n);

        gn = igamma(1-3/n,(xn/De).^n );

        e1 = exp(-(xmin/De)^n);
        e2 = exp(-(xmax/De)^n);      
        g1 = igamma(1-3/n,(xmin/De)^n);
        g2 = igamma(1-3/n,(xmax/De)^n );

        %g = real(g);
        g1 = real(g1);
        g2 = real(g2);

        A = ( e1-e2)/( g1-g2 );

        Fm(:) = (e1 - em)/(e1 - e2);
        fm(:) =  n / De^n * xm.^(n-1) .* em /(e1-e2);
        %fn(:) = n / De^n * xn.^(n-4) .* en /(e1-e2) * A;
        fn(:) = n / De^(n-3) * xn.^(n-4) .* en /(e1-e2) * A;
        Fn(:) = (g1-gn)/(g1-g2);

        % PD
    case 'PD'
        for i = 1:length(xm)
            if xm(i)<=xmin
                Fm(i) = 0;
            end

            if xm(i)>xmin
                Fm(i) = 1-(xmin./xm(i)).^n;
            end

            if xn(i)<=xmin
                Fn(i) = 0;
            end

            if xn(i)>xmin
                Fn(i) = 1-xn(i).^(-n-3)/xmin^(-n-3);
            end
        end

        fm(1) = Fm(1)-0;
        fn(1) = Fn(1)-0;
        for i = 2:length(xm)
            fm(i) = Fm(i) - Fm(i-1);
            fn(i) = Fn(i) - Fn(i-1);
        end

        % GSD
    case 'GSD'
        for i = 1:length(xm)
            if xm(i)>xmax
                Fm(i) = 1;
            end
            if xm(i)<=xmax
                Fm(i) = (xm(i)/xmax).^n;
            end

            if xn(i)>xmax
                Fn(i) = 1;
            end
            if xn(i)<=xmax
                Fn(i) = (xn(i)/xmax).^(n-3);
            end

        end

        fm(1) = Fm(1)-0;
        fn(1) = Fn(1)-0;
        for i = 2:length(xm)
            fm(i) = Fm(i) - Fm(i-1);
            fn(i) = Fn(i) - Fn(i-1);
        end

        % RRD
    case 'RRD'
        ee = exp(-(xm/De).^n);
        % A = k^3/igamma(1-3/n,0);
        % A = real(g);

        Fm(:) = 1-ee;
        % f(:,2) =  n / k^n * x.^(n-1) .* ee;
        % f(:,3) = n / k^n * x.^(n-4) .* ee * A;
        Fn(:) = 1-igamma(1-3/n,(xn/De).^n)/igamma(1-3/n,0);

        fm(1) = Fm(1)-0;
        fn(1) = Fn(1)-0;
        for i = 2:length(x)
            fm(i) = Fm(i) - Fm(i-1);
            fn(i) = Fn(i) - Fn(i-1);
        end

end
PSD_the = [Fm fm Fn fn];

end