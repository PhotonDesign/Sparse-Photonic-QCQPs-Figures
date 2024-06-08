function name = tName(varStr,varVal,headStr)
    N = length(varStr);
    name = headStr;
    for i = 1:N
        name_add = [varStr{i}, ' = ', num2str(varVal{i},4)];
        name = [name name_add ', '];
    end
    name = name(1:end-2);
end