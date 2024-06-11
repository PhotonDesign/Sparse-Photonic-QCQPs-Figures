function T = t2BTTB(t)
    % Create a block Toeplitz matrix with Toeplitz blocks (BCCB) T of
    % dimension nm by nm:
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

    % dimensions
    [N, M] = size(t);
    n = (N + 1) / 2;
    m = (M + 1) / 2;
    
    % create 2m - 1 Toeplitz blocks
    Ti = cell(1,M);
    for i = 1:M
        Ti{i} = t2Toeplitz(t(:,i));
    end

    % create the nm by nm BTTB matrix
    T = zeros(n*m, n*m);
    for i = 1:m % row index in block-level
        idx_row = ((i-1)*n + 1) : (i*n); % row index in element-level
        for j = 1:m % column index in block-level
            idx_col = ((j-1)*n + 1) : (j*n); % column index in element-level
            T(idx_row, idx_col) = Ti{i-j+m};
        end
    end
end


