function[D] = get_all_D(N)


D = cell(2*N,1); % ALL D matrices, real and complex

Dtmp = eye(N);

count = 1;

for kk = 1:2:length(D)
    
    D{kk} = diag(Dtmp(count,:));
    D{kk+1} = 1i*D{kk};
    
    count = count+1;
    
end