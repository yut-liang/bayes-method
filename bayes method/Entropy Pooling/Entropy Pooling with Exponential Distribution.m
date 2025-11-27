addpath('C:\Aorda\PSG\lib')
savepath 

clc
clear

% Prior distribution
rng(1)
% Analytical representation
N=2; % market dimension
lambda=1.5;  %Exponential Distribution
Mu=(1/lambda)*ones(N,1);
Sigma=diag(1/lambda^2 * ones(N,1));

% Numerical representation (panel)
J=5000;
p = ones(J,1)/J;

dd = exprnd(1/lambda, J/2, N);  
X = exprnd(1/lambda, J, N);  

% Investor views
%1 X = (1/2)Y
Q=[1 -2];
Mu_Q=0;
G = Q;
Sigma_G=.5^2

%2 Var(Y) = 4 Var(X)
SecMom = [1 -4] * Sigma * [1; -4] + Sigma_G;

% Posterior distribution

% Analytical posterior as the solution of relative entropy minimization
Mu_=Mu+Sigma*Q'*inv(Q*Sigma*Q')*(Mu_Q-Q*Mu);
Sigma_=Sigma+(Sigma*G')*( inv(G*Sigma*G')*Sigma_G*inv(G*Sigma*G') -inv(G*Sigma*G') )*(G*Sigma);

% Black-Litterman posterior
Mu_bl = inv((inv(Sigma)+Q'*inv(Sigma_G)*Q))*(inv(Sigma)*Mu + Q'*inv(Sigma_G)*Mu_Q);
Sigma_bl = inv(inv(Sigma)+Q'*inv(Sigma_G)*Q);

% Numerical posterior distribution of entropy pooling
Aeq = ones(1,J);  % constrain probabilities to sum to one...
beq=1;
QX = X*Q';
Aeq = [Aeq   % ...constrain the first moments...
    QX'];
beq=[beq
    Mu_Q];

GX = X*G';
for k=1:size(G,1)
    for l=k:size(G,1)
        Aeq = [Aeq
            (GX(:,k).*GX(:,l))'];
        beq=[beq
            SecMom(k,l)];
    end
end

tic
% Entropy pooling using riskprog
H=p';
lb = ones(J,1)*(10^(-13));
stroptions.Linearization = 'On';
stroptions.Solver = 'BULDOZER';
[p_psg, fval, status, output] = riskprog('entropyr', [], H, [], [], [], [], [], Aeq, beq, lb, [], [], stroptions);
toc


% Figures

[J,N]=size(X);
NBins=round(10*log(J));
for n=1:N
    figure 
    
    % Set ranges
    xl=min(X(:,n));
    xh=max(X(:,n));
    x=[xl : (xh-xl)/100 : xh];
    
    % Posterior entropy pooling
    pHist(X(:,n),p_psg,NBins);
   

    % Posterior analytical
    hold on
    lambda_post = 1 / Mu_(n); 
    plot(x, exppdf(x, lambda_post), 'linewidth', 2, 'color', 'r');

    % Black-Litterman posterior 
    hold on
    lambda_bl = 1 / Mu_bl(n); 
    plot(x, exppdf(x, lambda_bl), 'linewidth', 2, 'color', 'g');

    % Prior analytical
    hold on
    plot(x, exppdf(x, lambda), 'linewidth', 2, 'color', 'b');
    
    
    xlim([xl xh])
    legend('entropy pooling', 'analytical', 'BL', 'prior')
end

