# Diamond_Simulation
Matlab simulation for model developed upon Diamond model. Thoughts and derivations can be found at https://github.com/YinshanS/Thoughts-on-Diamond 


## Contents:

### [ql_u_2s_rf.m](https://github.com/YinshanS/Diamond_Simulation/blob/master/ql_u_2s_rf.m)
calculates risk-free interest rate corresponding to different
debt level, assuming
1. constant young marginal utility
2. Cobb-douglas production
3. 2-state shock 

### [ql_u_2s.m](https://github.com/YinshanS/Diamond_Simulation/blob/master/ql_u_2s.m)
assumes 
1. linear young utility and log old utility,
    U=cy+beta*ln(co)+p(public goods)
2. cobb-douglas production
    Y_t=A(1+epsilon_t)K_{t-1}^alpha
3.  2-state uncertainty
    epsilon_t=+-sigma w/ prob 1/2
4. Fixed debt amount(if enough for paying back)

and calculates the debt amount and welfare effects of the debt policy.

### [log_u_2s.m](https://github.com/YinshanS/Diamond_Simulation/blob/master/log_u_2s.m)
assumes 
1. log utility,
    U=(1-beta)*ln(cy)+beta*ln(co)+p(public goods)
2. cobb-douglas production
    Y_t=A(1+epsilon_t)K_{t-1}^alpha
3.  2-state uncertainty
    epsilon_t=+-sigma w/ prob 1/2
4. Fixed debt amount(if enough for paying back)

and calculates welfare effect and debt amount of the proposed fiscal policy.
