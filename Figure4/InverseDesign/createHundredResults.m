
for i=11:40
    [mHist,xHist,fc,x,fval] = test3c();
    save(['HundredResultsGradDescent/data',num2str(i),'.mat'],'fc','mHist','xHist');
end