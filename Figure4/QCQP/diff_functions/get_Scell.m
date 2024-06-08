function Scell = get_Scell(S,F,nc,bags)
    Scell = cell(length(F));
    for i = 1:(2*nc)
        Scell{i} = S(i,i);
    end
    start_ind = 2*nc+1;
    for i = 1:length(bags)
        end_ind = length(bags{i})-1+start_ind; 
        %Blocks are same size as bags
        %Size of bags - 1 but + 1 because of extended submat
        Scell{i} = S(start_ind:end_ind,start_ind:end_ind);
        start_ind = end_ind + 1;
    end
    Scell{end} = S(start_ind:end,start_ind:end);
end