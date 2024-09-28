
for i=1:50
    [mHist,xHist,fc,xval,fval] = inverse_designMMA();
    save(['MMAFiftyResultsRandStart/data',num2str(i),'.mat'],'fc','mHist','xHist');
end