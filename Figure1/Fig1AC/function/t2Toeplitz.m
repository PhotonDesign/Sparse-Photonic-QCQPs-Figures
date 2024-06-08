function T = t2Toeplitz(t)
    % input: a 2n-1 by 1 vector 
    %
    % t = [ t_{-(n-1)} ]
    %     [     :      ]
    %     [   t_{0}    ]
    %     [   t_{1}    ]
    %     [     :      ]
    %     [  t_{n-1}   ] 
    %
    % output: a n by n Toeplitz matrix
    %
    % T = [ t_{0}   t_{-1}  ...  t_{-(n-1)} ]
    %     [ t_{1}   t_{0}   ...  t_{-(n-2)} ]
    %     [ ...     ...     ...     ...     ]
    %     [ t_{n-1} t_{n-2} ...    t_{0}    ]
  
    Nt = length(t);
    n = (Nt + 1) / 2;
    c = t(n : 2*n-1);
    r = t(n : -1 : 1);
    T = toeplitz(c,r);
end