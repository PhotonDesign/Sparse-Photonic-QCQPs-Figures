function [Qext, Qabs, Qsca, r, t] = Qana_Q_film(eps, th1, h)
    
    n = sqrt(eps);

    n1 = 1;
    n2 = n;
    n3 = 1;

    % calculation
    th2 = asin(n1/n2*sin(th1)); % radius
    th3 = asin(n1/n3*sin(th1)); % radius

    b = n2*h*cos(th2);

    p1 = n1 * cos(th1);
    p2 = n2 * cos(th2);
    p3 = n3 * cos(th3);

    r12 = (p1-p2) ./ (p1+p2);
    t12 = 2*p1 ./ (p1+p2);
    r23 = (p2-p3) ./ (p2+p3);
    t23 = 2*p2 ./ (p2+p3);

    r = (r12 + r23.*exp(1j*2*b)) ./ (1+r12.*r23.*exp(1j*2*b));
    t = (t12.*t23.*exp(1j*b)) ./ (1+r12.*r23.*exp(1j*2*b));

    t2 = t ./ exp(1j*h*cos(th1)); % normalize by incident phase
    
    
    
    Qext = 2*real(1-t2)*cos(th1);
    Qabs = (1-abs(t2).^2-abs(r).^2)*cos(th1);
    Qsca = (abs(r).^2+abs(t2-1).^2)*cos(th1);
end