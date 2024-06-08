function [St, bt, ct] = c2r_QCQP(S, b, c)
    % convert a real-valued QCQP with complex-valued variable x
    %   x'*S*x + real(b'*x) + c
    % to a real-valued QCQP with real-valued variable xt = [Re(x), Im(x)]
    %   xt'*St*xt + real(bt'*xt) + ct
    St = [S 1j*S; -1j*S S];
    bt = [b; -1j*b];
    ct = c;
end
