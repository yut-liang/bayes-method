%clc; 
%clear;

% 1-5
%Load data:
data = readtable('C:\Users\ASUS\Desktop\AMS_ Homework8 _data.csv');

prices = table2array(data(:,2:5));  

% log return
returns = diff(log(prices));  


% regression
H{1} = returns(:,2:4);  % MSFT, GOOGL, AMZN
c{1} = returns(:,1);    % AAPL


count = 1;
iargstruc_arr(count) = matrix_pack('matrix_s',H{1},[],c{1});

%Define problem statement for CVaR regression:
problem_statement_cvar_2 = sprintf('%s\n',...
'minimize',...
'  cvar2_err(0.8,matrix_s)',...
'Value:',...
'  L(matrix_s)',...
' ');

%Optimize problem using mpsg_solver function:
[solution_str_2, outargstruc_arr_2] = mpsg_solver(problem_statement_cvar_2, iargstruc_arr);

tail_los_full_cvar = outargstruc_arr_2(3).values(:,2); 
intercept_cvar = outargstruc_arr_2(1).values(1);  

tail_los_cvar = tail_los_full_cvar + intercept_cvar;

cvar_0_80 = functionvalue('cvar_risk_g', 0.8, tail_los_cvar, [], [], 1);

rep_text = 'CVaR regression with parameter alpha = 0.8';
rep_text = sprintf('%s \n intercept = %f, CVaR = %f',rep_text,intercept_cvar,cvar_0_80);


var_080 = functionvalue('var_risk_g', 0.8, tail_los_cvar, [], [], 1);
loss_null_var = tail_los_cvar - var_080;
loss_positive = loss_null_var(loss_null_var>0);

%Plot histogram of tail of loss of the CVaR regression
h2 = figure;
histogram(loss_positive,30,'FaceColor',[0.5,0.5,0.5])

%Estimate parameters of GPD
params_cvar = tsallis_harmonik_params_loss(loss_positive);

%Report for parameters estimation of GPD for tail of loss of the CVaR
%regression
rep_text = sprintf('%s\n Mu = %f\n',rep_text,params_cvar(1).mua);
rep_text = sprintf('%s Kappa MLE = %f\n',rep_text,params_cvar(1).kappa);
rep_text = sprintf('%s Kappa Harmonic = %f',rep_text,params_cvar(2).kappa);

disp('    ')
disp(rep_text)





% 6
% choose day 100 return
day_idx = 100;

% find MSFT, GOOGL, AMZN return
x_new = returns(day_idx, 2:4);   % 1x3 

% regression
full_beta = outargstruc_arr_2(1).values;   % e.g. [intercept, b1, b2, b3]
intercept = full_beta(1);
beta = full_beta(1:end);         % beta: 3x1
beta = beta(:);                

predicted_loss = x_new * beta + intercept;

fprintf('Predicted loss of AAPL on day %d = %.6f\n', day_idx, predicted_loss);






% 7
num_days = size(returns,1);   
predicted_losses = zeros(num_days,1); 

% regression parameter
full_beta = outargstruc_arr_2(1).values;
intercept = full_beta(1);
beta = full_beta(1:end);
beta = beta(:);

% everyday predict loss
for t = 1:num_days
    x_row = returns(t, 2:4);     % MSFT, GOOGL, AMZN return
    predicted_losses(t) = x_row * beta + intercept;
end

%  VaR and Tail Loss 
alpha = 0.8;
VaR_pred = quantile(predicted_losses, 1 - alpha);  % 20% quantile
tail_losses_pred = predicted_losses(predicted_losses <= VaR_pred);

% graph
figure;
histogram(tail_losses_pred, 30, 'FaceColor', [0.4 0.4 0.8])
title(sprintf('Tail Loss Distribution (Predicted, alpha = %.2f)', alpha))
xlabel('Predicted Loss')
ylabel('Frequency')


