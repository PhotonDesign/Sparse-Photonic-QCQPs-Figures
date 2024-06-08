

function [y] = sqrt_k(x)
    a = angle(x);
    a(a<0) = a(a<0) + 2*pi;
    y = sqrt(abs(x)) .* exp(1j*a/2);

end