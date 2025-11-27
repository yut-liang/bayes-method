% AMZN, NVDA, TSLA
mu = [0.001594; 0.003206; 0.000728]; 
cov_matrix = [0.000334, 0.000241, 0.000230;
              0.000241, 0.001033, 0.000316;
              0.000230, 0.000316, 0.001278];  

num_assets = length(mu); 
N = 250;

% from python: market portfolio's μ and σ
mu1 = 0.000786;
cov1 = 0.000063;

% find sharp
mean_spy = mu1 * 250;  % SPY annual return
std_spy = sqrt(cov1) * sqrt(250); 
shpr = mean_spy / std_spy;

wtsMarket = [0.344533; 0.193961; 0.461505];  

% INVESTOR VIEWS
n_views = 3;  
P_matrix = zeros(n_views, num_assets);
q = zeros(n_views, 1);
omega = zeros(n_views);

% INVESTOR VIEWS
n_views = 3;  
P_matrix = zeros(n_views, num_assets);
q = zeros(n_views, 1);
omega = zeros(n_views);

% #1 Absolute view: Expected annualized return of AMZN is greater than 0.005
P_matrix(1, 1) = 1; 
q(1) = 0.005;
omega(1, 1) = 0.001;

% #2 Relative view: NVDA has an annualized return 0.002 higher than TSLA
P_matrix(2, 2) = 1;  % NVDA
P_matrix(2, 3) = -1; % TSLA 
q(2) = 0.002;  
omega(2, 2) = 0.02;

% #3 Absolute view: Expected annualized return of TSLA is greater than 0.0012
P_matrix(3, 3) = 1;  % TSLA
q(3) = 0.0012;  
omega(3, 3) = 0.05;

conv_coef = 1 / 250; 
q = q * conv_coef;
omega = omega * conv_coef;


P_eq = delta * cov_matrix * mu;  

disp('P_matrix:');
disp(P_matrix);
disp('q (daily returns):');
disp(q);
disp('omega:');
disp(omega);

%Solve portfolio optimization problem
numAssets = num_assets;
LB = zeros(1, numAssets);
Aeq = ones(1, numAssets);
Beq = 1;

% Define uncertainty tau
tau = [0.3, 0.5, 0.9];

% Find delta (只保留一行)
delta = shpr / sqrt(wtsMarket' * cov_matrix * wtsMarket);

% 计算 Black-Litterman 调整后的均值和协方差
mean_mu_bl = zeros(num_assets, length(tau));
cov_mu_bl = cell(1, length(tau));

for i = 1:length(tau)
    mean_mu_bl(:, i) = inv(inv(tau(i) * cov_matrix) + P_matrix' * inv(omega) * P_matrix) * ...
        (inv(tau(i) * cov_matrix) * P_eq + P_matrix' * inv(omega) * q);
    cov_mu_bl{i} = inv(inv(tau(i) * cov_matrix) + P_matrix' * inv(omega) * P_matrix);
end

disp('Black-Litterman Adjusted Expected Returns (Annualized):');
for i = 1:length(tau)
    fprintf('Tau = %.1f\n', tau(i));
    disp(mean_mu_bl(:, i) * 250);
end

% Optimal portfolio allocation
wts_alloc = zeros(num_assets, length(tau));
for i = 1:length(tau)
    wts_alloc(:, i) = inv(cov_matrix + cov_mu_bl{i}) * mean_mu_bl(:, i) / delta;
end

figure;
hold on;
bar(wts_alloc, 'grouped');
legend('Tau = 0.3', 'Tau = 0.5', 'Tau = 0.9');
title('Optimal Portfolio Weights with Black-Litterman');
xlabel('Assets');
ylabel('Weights');
xticks(1:num_assets);
xticklabels({'AMZN', 'NVDA', 'TSLA'}); 
hold off;









%Estimate efficient portfolio that maximizes Sharpe ratio
port = Portfolio('NumAssets', num_assets, 'lb', 0, 'budget', 1, 'Name', 'Mean Variance');
port = setAssetMoments(port, mu, cov_matrix);
wts = estimateMaxSharpeRatio(port);

wtsBL = zeros(num_assets, length(tau));
for i = 1:length(tau)
    portBL = Portfolio('NumAssets', num_assets, 'lb', 0, 'budget', 1, ...
        'Name', ['Mean Variance with Black-Litterman (Tau = ', num2str(tau(i)), ')']);
    portBL = setAssetMoments(portBL, mean_mu_bl(:, i), cov_matrix + cov_mu_bl{i});
    wtsBL(:, i) = estimateMaxSharpeRatio(portBL);
end


asset_names = {'AMZN', 'NVDA','TSLA'};
figure;

subplot(1, length(tau) + 1, 1); 
idx_mv = wts > 0.001;
pie(wts(idx_mv), asset_names(idx_mv));
title('Mean-Variance Portfolio Weights');

for i = 1:length(tau)
    subplot(1, length(tau) + 1, i + 1); 
    idx_BL = wtsBL(:, i) > 0.001;
    pie(wtsBL(idx_BL, i), asset_names(idx_BL));
    title(['BL Portfolio Weights (Tau = ', num2str(tau(i)), ')']);
end

result_optim_weights = table(asset_names', wts, wtsBL(:, 1), wtsBL(:, 2), wtsBL(:, 3), 'VariableNames', ...
    {'Asset_Name', 'Mean_Variance', 'BL_Tau_0.3', 'BL_Tau_0.5', 'BL_Tau_0.9'});
disp(result_optim_weights);

