function [pmlx, pmly, pmlz] = PML(Dim,h,BC)
    N = round(Dim(1)/h);
    M = round(Dim(2)/h);
    X = 0:h:Dim(1);
    Y = 0:h:Dim(2); 
    Xs = h/2:h:Dim(1); %staggered x grid for matx
    Ys = h/2:h:Dim(2); %staggered y grid for maty
    
    %take only N, M grids
    X = X(1:N);
    Y = Y(1:M).';
    Xs = Xs(1:N);
    Ys = Ys(1:M).';
    
    switch BC{1}{1}
        case 'pml'
            sig_x = @(x) BC{1}{2}(2)*(x-Dim(1)+BC{1}{2}(1)).*(x>=(Dim(1)-BC{1}{2}(1))) - BC{1}{2}(2)*(x-BC{1}{2}(1)).*(x<=BC{1}{2}(1));
        case 'periodic'
            %no pml
            sig_x = @(x) zeros(1,N);
        otherwise
            warning('invalid x boundary condition');
    end
            
    switch BC{2}{1}
        case 'pml'
            sig_y = @(y) BC{2}{2}(2)*(y-Dim(2)+BC{2}{2}(1)).*(y>=(Dim(2)-BC{2}{2}(1))) - BC{2}{2}(2)*(y-BC{2}{2}(1)).*(y<=BC{2}{2}(1));
        case 'periodic'
            %no pml
            sig_y = @(y) zeros(M,1);
        otherwise
            warning('invalid x boundary condition');
    end
    
    
    pmlx = ((1 + 1i*sig_y(Y))*(1./(1 + 1i*sig_x(Xs))));
    pmly = (1./(1 + 1i*sig_y(Ys))*((1 + 1i*sig_x(X))));
    pmlz = ((1 + 1i*sig_y(Y))*(1 + 1i*sig_x(X)));
end
