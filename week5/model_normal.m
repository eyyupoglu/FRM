function [VAR99, VAR99_EWMA] = model_normal(v, alpha, w)
N = (size(v,1)-1);

HistoricalReturns = diff(v)./v(1:end-1,:);

CorrMatrix = corr(HistoricalReturns); %Table 14.3
CovMatrix = cov(HistoricalReturns); %Table 14.4 (slight deviations on 6th decimal)

% NOTE: Slightly different values due to (probably) rounding errors in
% Matlab/Excel or use of difference in variance/covariance formulas

PF_variance = alpha'*CovMatrix*alpha; %equation 14.3, measured in $000s
returns_variance = w'*CovMatrix*w;
VAR99 = norminv(0.99)*sqrt(PF_variance);
VAR99_returns = norminv(0.99)*sqrt(returns_variance);
ES99 = sqrt(PF_variance)*exp(-norminv(0.99)^2/2)/sqrt(2*pi)/0.01;
VAR99 - VAR99_returns*10000 % deviation on 15th decimal or so

VAR95 = norminv(0.95)*sqrt(PF_variance);
VAR95_returns = norminv(0.95)*sqrt(returns_variance);
ES95 = sqrt(PF_variance)*exp(-norminv(0.95)^2/2)/sqrt(2*pi)/0.05;
VAR95 - VAR95_returns*10000 % deviation on 15th decimal or so






lambda = 0.94;
VariancesEWMA = CovMatrix; %starting point - it does not matter much when using a long time series. Try to start in Identity matrix instead

for i = 2:N+1

   VariancesEWMA(:,:,i) = lambda*VariancesEWMA(:,:,i-1) + (1-lambda)*HistoricalReturns(i-1,:)'*HistoricalReturns(i-1,:);
    
end

% Table 14.6
FinalVols(1,:) = 100*diag(CovMatrix).^0.5;
FinalVols(2,:) = 100*diag(VariancesEWMA(:,:,end)).^0.5;

SD = diag(VariancesEWMA(:,:,end)).^0.5;

% Table 14.7
CorrMatrixEnd = diag(1./SD)*VariancesEWMA(:,:,end)*diag((1./SD));

PF_variance_EWMA = alpha'*VariancesEWMA(:,:,end)*alpha;
VAR99_EWMA = norminv(0.99)*sqrt(PF_variance_EWMA);
ES99_EWMA = sqrt(PF_variance_EWMA)*exp(-norminv(0.99)^2/2)/sqrt(2*pi)/0.01;

VAR95_EWMA = norminv(0.95)*sqrt(PF_variance_EWMA);
ES95_EWMA = sqrt(PF_variance_EWMA)*exp(-norminv(0.95)^2/2)/sqrt(2*pi)/0.05;


end