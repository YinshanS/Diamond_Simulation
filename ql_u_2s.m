%% Notes
%{
This .m file assumes 
1. linear young utility and log old utility,
    U=cy+beta*ln(co)+p(public goods)
2. cobb-douglas production
    Y_t=A(1+epsilon_t)K_{t-1}^alpha
3.  2-state uncertainty
    epsilon_t=+-sigma w/ prob 1/2
4. Fixed debt amount(if enough for paying back)

and calculates the debt amount and welfare effects of the debt policy.
%}

%% Prep
clear
close all
clc
%% Parameters
N = 10000;
T = 200;
alpha = 0.35;
beta = 0.495;
K0 = beta;
A0 = 10/(3*beta^(alpha-1));%A in production function, set to make Y/K=10/3
V=10;
Q=8;
%% Simulation
dU = zeros(N,Q+1,V);
blowup = zeros(Q+1,V);
for s=1:N
    shock = binornd(1,0.5,T,1)*2-1;
    %A = ((binornd(1,0.5,T,1)*2-1)*sigma+1)*A0; %2-state shock
for v=1:V
    D = 0.001*v;
for q=1:Q+1
    sigma = 0.4+0.05*(q-1);
    A = (shock*sigma+1)*A0;
    [dU(s,q,v),blowup(q,v)] = simu(T,sigma,alpha,beta,K0,A0,A,D,blowup(q,v));
end
end
end
save('ql_u_2s.mat');
dU1=squeeze(mean(dU,1));
b=bar3(dU1);
title({'Welfare Effect of Constant Debt Policy with Tax and Subsidy','(Constant Young Marginal Utility, 2-state Uncertainty)','(\alpha=0.35, \beta=0.495, T=200, A=10\beta^{1-\alpha}/3)'});
ylabel('Standard Deviation of 2-state Shock','Position',[-1.86324020745124,3.6469991034619,-0.002477681111464]); 
xlabel('Debt Amount') ;
xticklabels({'0.001','0.002','0.003','0.004','0.005','0.006','0.007','0.008','0.009','0.01'});
yticklabels({'0.4','0.45','0.5','0.55','0.6','0.65','0.7','0.75','0.8'});
zlabel('Welfare Change');
savefig('Welfare Effect Constant Debt.fig');

%blowup!



%% Function
function[dU,blowup] = simu(T,sigma,alpha,beta,K0,A0,A,D,blowup)
%without debt:
Y = zeros(T,1); 
W = zeros(T,1);
U = zeros(T-1,1);
P = zeros(T,1);%government purchased public good
K = zeros(T,1);
r = zeros(T,1);%risky return rate
cy = zeros(T,1);
co = zeros(T,1);

Y(1) = A(1)*K0^alpha;
W(1) = (1-alpha)*Y(1);
K(1) = min(beta,W(1));%invest beta or all wage(if wage<beta)
cy(1) = W(1)-K(1);
co(1) = Y(1)-W(1);
r(1) = alpha*A(1)*K0^(alpha-1);

%with debt:
% variable1 represents variable with debt
Y1 = zeros(T,1); 
W1 = zeros(T,1);
U1 = zeros(T-1,1);
P1 = zeros(T,1);%government purchased public good
K1 = zeros(T,1);
r1 = zeros(T,1);%risky return rate
cy1 = zeros(T,1);
co1 = zeros(T,1);

D1 = zeros(T,1);%debt
D1(1) = D;
rsb1 = zeros(T,1);%risk-free interest rate
tau = zeros(T,1);

Y1(1) = A(1)*K0^alpha;
W1(1) = (1-alpha)*Y1(1);
K1(1) = min(beta-D,W1(1)-D);%invest beta or all wage(if wage<beta)
cy1(1) = W1(1)-K1(1)-D;
co1(1) = Y1(1)-W1(1)+D;
r1(1) = A(1)*K0^(alpha-1);
tau(1) = r1(1)/r(1)-1;
%at time 1 the government gives all debt income as tansfer to the old
%in subsequent periods the government repays debt first and use surplus to
%provide public goods to the old

for t=2:T
    Y(t) = A(t)*K(t-1)^alpha;
    W(t) = (1-alpha)*Y(t);
    K(t) = min(beta,W(t));
    cy(t) = W(t)-K(t);
    co(t) = alpha*Y(t);
    U(t-1) = cy(t-1)+beta*log(co(t))+P(t);
    r(t) = alpha*Y(t)/K(t-1);  
    
    Y1(t) = A(t)*K1(t-1)^alpha;
    r1(t) = alpha*Y1(t)/K1(t-1);
    tau(t) = 1-r(t)/r1(t);
    rsb1(t) = (1-tau(t))*(2*alpha*A0*(beta-D1(t-1))^alpha*(1-sigma^2))/(beta-2*D1(t-1)+sqrt((beta-2*D1(t-1))^2+4*(1-sigma^2)*D1(t-1)*(beta-D1(t-1))));
    D1(t)=max(D, D1(t-1)*rsb1(t));%constant debt
    %D1(t)=D1(t-1)*rsb1(t)+W(t-1)-W1(t-1);
    %D1(t)=D1(t-1)*rsb1(t)+(W(t-1)-W1(t-1))*r(t);
    W1(t) = Y1(t)-r(t)*K1(t-1);
    %W1(t) = (1-alpha)*Y1(t);
    K1(t) = min(beta-D1(t),W1(t)-D1(t));
    if isreal(rsb1(t))== 0 || K1(t)<=0
       disp(['negative capital at'... 
        num2str(t)]);
        blowup = blowup+1;
       break;
    end
    cy1(t) = W1(t)-K1(t)-D1(t);
    %co1(t) = K1(t-1)*r(t)+D1(t);
    co1(t) = K1(t-1)*r(t)+D1(t-1)*rsb1(t);
    %co1(t) = Y1(t)-W1(t)+D1(t-1)*rsb1(t);
    P1(t) = max(0,D1(t)-D1(t-1)*rsb1(t));
    U1(t-1) = cy1(t-1)+beta*log(co1(t))+P1(t);
end
dU = mean(U1-U);
end
%{
disp(['Y='...
    num2str(Y')]);
disp(['W='...
    num2str(W')]);
disp(['K='...
    num2str(K')]);
disp(['cy='...
    num2str(cy')]);
disp(['co='...
    num2str(co')]);
disp(['r='...
    num2str(r')]);
disp(['U='...
    num2str(U')]);

disp(['Y1='...
    num2str(Y1')]);
disp(['W1='...
    num2str(W1')]);
disp(['K1='...
    num2str(K1')]);
disp(['rsb1='...
    num2str(rsb1')]);
disp(['cy1='...
    num2str(cy1')]);
disp(['co1='...
    num2str(co1')]);
disp(['D1='...
    num2str(D1')]);
disp(['P1='...
    num2str(P1')]);
disp(['U1='...
    num2str(U1')]);
disp(['U1-U='...
    num2str(U1'-U')]);
disp(['blowup='...
    num2str(blowup)]);

%It is normal that for some generation the utility might be decreased.
%However they are maximizing expected utility. The realized value might be
%undesirable sometimes, but overall they achieve an improvement in expected utility.
%}

