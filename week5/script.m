% HullData = HullDataExtended.Data;
% v = [HullData(1:501,1) HullData(1:501,2).*HullDataExtended.HullFX(:,1) HullData(1:501,3:4)./HullDataExtended.HullFX(:,2:3)]; % Table 13.2


% [VAR99, VAR99_EWMA] = model_normal(v);
clc
clear
load 'StockFXData.mat'


alpha = [2000 2000 2000 2000 2000]';
w = ([2000/10000 2000/10000 2000/10000 2000/10000 2000/10000])';
N = 250;

v = FXCopy.Data(:, 1:5) .* FXCopy.USEU;

% v = [HullDataExtended.USData(end-500-i:end-i,:)]; % Table 13.2
[VAR99_1, VAR99_EWMA] = model_normal(v, alpha, w);



%% 1.2


levels = [FXCopy.USEU FXCopy.Data(:, 1:5)];
r_f = (levels(2:end, :) - levels(1:end-1, :)) ./ levels(1:end-1, :);

covariance = cov(r_f);
% w_map = ([10000/10000 2000/10000 2000/10000 2000/10000 2000/10000 2000/10000])';

w_map = [10000 2000 2000 2000 2000 2000]';
PF_variance = w_map' * covariance * w_map;
VAR99_2 = norminv(0.99)*sqrt(PF_variance);

%% 1.3
levels = [FXCopy.USEU FXCopy.Data];

r_f = (levels(2:end, 2:end-1) - levels(1:end-1, 2:end-1)) ./ levels(1:end-1, 2:end-1);

r_f_i = (levels(2:end, end) - levels(1:end-1, end)) ./ levels(1:end-1, end);

beta_1 = fitlm(r_f_i, r_f(:, 1)).Coefficients.Estimate(2);
beta_2 = fitlm(r_f_i, r_f(:, 2)).Coefficients.Estimate(2);
beta_3 = fitlm(r_f_i, r_f(:, 3)).Coefficients.Estimate(2);
beta_4 = fitlm(r_f_i, r_f(:, 4)).Coefficients.Estimate(2);
beta_5 = fitlm(r_f_i, r_f(:, 5)).Coefficients.Estimate(2);

% (r_f * r_f')^(-1) * r_f .* r_f_i

beta = [beta_1 beta_2 beta_3 beta_4 beta_5];
w = w_map(2:end).' * beta.';

PF_variance = w' * cov(r_f_i) * w;
VAR99_3 = norminv(0.99)*sqrt(PF_variance);

%% 1.4

levels = [FXCopy.USEU FXCopy.Data(:, :) .* FXCopy.USEU];

r_f = (levels(2:end, 2:end-1) - levels(1:end-1, 2:end-1)) ./ levels(1:end-1, 2:end-1);
r_f_i = (levels(2:end, end) - levels(1:end-1, end)) ./ levels(1:end-1, end);

w = ([2000/10000 2000/10000 2000/10000 2000/10000 2000/10000])';
r_portf = r_f * w;

beta_1 = fitlm(r_portf, r_f_i).Coefficients.Estimate(2);


beta = [beta_1];
w = 10000 * beta.';
PF_variance = w * cov(r_portf) * w;
VAR99_4 = norminv(0.99)*sqrt(PF_variance);




