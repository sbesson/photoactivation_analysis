
conditions(1).series(1).data
i=1;
iEvent = conditions(1).series(1).events(2).index(i);
y = conditions(1).series(1).data(iEvent:end,i);
x = conditions(1).series(1).times(iEvent:end);


fitFun = @(b,x)(b(1) + b(2) .* exp(- x ./ b(3)));     %Double-exponential function for fitting
bInit = [y(end) y(1)-y(end) 10];
%Fit function to ratio timeseries
fitOptions = statset('Robust','on','MaxIter',500,'Display','off');
[bFit,resFit,jacFit,covFit,mseFit] = nlinfit(x(:) - x(1),...
    y(:),fitFun,bInit,fitOptions);
%Get confidence intervals of fit and fit values
[fitValues,deltaFit] = nlpredci(fitFun,x(:),bFit,resFit,'covar',covFit,'mse',mseFit);

%Ch