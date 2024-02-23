% ODE4 Algorithm from:
% John Keevil (2024). ODE4 gives more accurate results than ODE45,
% ODE23, ODE23s (https://www.mathworks.com/matlabcentral/fileexchange/59044-ode4-gives-more-accurate-results-than-ode45-ode23-ode23s),
% MATLAB Central File Exchange. Retrieved February 21, 2024.
clearvars
close all
clc
name_model = 'model'; 
parameters_model;
param_model = [bP bI bA th dE dP sig nI nA z gI gA gH gQ aI aH Hbed Nd Npop X0];
parameters_cost;
param_cost = [w r VSL Xd p_o p_y cH ts cV cU];
clk = clock;
title_format = 'day_%d-%d-%d__hour_%d-%d__model_%s.mat';
title_file = sprintf(title_format,clk(1),clk(2),clk(3),clk(4),clk(5),name_model);
clear clk title_format 
ts = [0 15 30 45 60 75 90 105 120 135 150];
Lzero = zeros(Nd+1,1);
Xfree = ode4(name_model, 0:Nd, X0, param_model, Lzero);
X0_table = [ts' Xfree((ts)+1,:)]; %Initial state conditions 
clearvars Lzero
Lambda = 0.00:0.02:1.00;
ub = 0.7*ones(1,length(time));
lb = zeros(1,length(time)); 
Nini = 1;
U0(:,Nini) = 0.7.*rand(366,1);
Nlk = size(U0,2);
parpool('local')
Nsim = size(Lambda,2);
Nts = size(ts,2);
Uopt = zeros(Nd+1,Nsim,Nts);
Xopt = zeros(Nd+1,length(X0),Nsim,Nts);
cost_val = zeros(Nsim,Nts);
exitflags = zeros(Nsim,Nts); 
for ii = 1:Nsim   
    Lam = Lambda(ii); %Lambda(:,ii);    
    for kk = 1:Nts 
        param_model(20:end) = Xfree(ts(kk)+1,:); 
        param_model(18) = Nd - ts(kk); 
        cost_functional = @(Lvect)cost_function(Lvect, param_model, param_cost, Lam, name_model);
        text = sprintf("Optimizing Lam combination %d/%d. Optimal control starts at day %d",ii,Nsim,ts(kk)); disp(text)
        U0_try = U0(ts(kk)+1:end);
        options = optimoptions(@fmincon,'Algorithm','sqp', 'Display', 'iter-detailed',...
    'UseParallel', true, 'MaxFunctionEvaluations',2000*length(U0_try));
        [Uopt_try,cost_val(ii,kk),exitflags(ii,kk)] = fmincon(cost_functional,U0_try,[],[],[],[],lb(ts(kk)+1:end),ub(ts(kk)+1:end),[],options);
        Uopt(:,ii,kk) = [zeros(ts(kk),1); Uopt_try];
        Xopt(:,:,ii,kk) = ode4(name_model, time, X0_table(1,2:end), param_model, Uopt(:,ii,kk));
        clc        
    end
end
save(title_file,'name_model','Lambda','Uopt','Xopt','U0','cost_val','exitflags',...
    'param_model','param_cost','X0_table'); 
delete(gcp) 
