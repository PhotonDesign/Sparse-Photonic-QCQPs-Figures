function y = BTTB_matvec(t, x)
    % Perform fast BTTB matrix-vector mutiplication y = Tx where T is a
    % nm by nm block toeplitz matrix with toeplitz blocks (BTTB) and x is
    % an nm by 1 vector.
    %
    % T = [ T_{0}    T_{-1}   ...   T_{-(m-1)}]
    %     [ T_{1}    T_{0}    ...   T_{-(m-2)}]
    %     [  ...     ...      ...      ...    ]
    %     [ T_{m-1}  T_{m-2}  ...     T_{0}   ] 
    %
    % Each block T_{i} is a Toeplitz matrix of dimension n by n,
    % which is completely defined by a 2n-1 by 1 vector t_{i} (see function
    % t2Toeplitz for detail.)
    %
    % The input t is a 2n-1 by 2m-1 matrix constructed by concatenating all the
    % t_{i} horizontally as such:
    %    
    % t = [     :        :        :        :        :    ]
    %     [     :        :        :        :        :    ]
    %     [ t_{-(m-1)}  ...      t_{0}    ...    t_{m-1} ]
    %     [     :        :        :        :        :    ]
    %     [     :        :        :        :        :    ]
    
    % note: x could also directly be a reshaped n by m matrix
    
    % dimensions
    [N, M] = size(t);
    n = (N + 1) / 2;
    m = (M + 1) / 2;
    
    % convert the nm by 1 input vector into a n by m matrix and add zero
    % paddings to stuff it to 2n-1 by 2m-1
    X = reshape(x, n, m);
    Xpp = zeros(N, M); % X with zero paddings
    Xpp(1:n, 1:m) = X;

    % construct a embedded circulant matrix by swapping rows and columns
    % of the input convolution matrix
    swap_row = [n:(2*n-1), 1:(n-1)];
    swap_col = [m:(2*m-1), 1:(m-1)];
    cpp = t(swap_row, swap_col); 

    % get output and strip down the paddings, reshape the resulting matrix 
    % to vector form 
    Ypp = ifft2(fft2(cpp) .* fft2(Xpp));
    Y = Ypp(1:n, 1:m);
    y = reshape(Y, n*m, 1);
end