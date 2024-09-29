function [Sr, br, cr] = c2r_con(S, b, c)
    ND = length(S);
    Sr = cell(1,2*ND);
    br = cell(1,2*ND);
    cr = cell(1,2*ND);
    for i = 1:ND
        Sr{2*i-1} = Re(S{i});
        Sr{2*i} = Im(S{i});
        br{2*i-1} = b{i};
        br{2*i} = 1j*b{i};
        cr{2*i-1} = real(c{i});
        cr{2*i} = imag(c{i});
    end
end
