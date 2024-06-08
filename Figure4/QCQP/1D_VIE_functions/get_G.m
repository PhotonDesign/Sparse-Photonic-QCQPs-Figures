function G = get_G(x1,x2,k,dx)
    G = (1i*k/2)*exp(1i*k*abs(x1-x2))*dx; %dx should be L/M
end