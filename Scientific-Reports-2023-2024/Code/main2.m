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
Lambda = 0:0.02:1; %51
ub = 0.7*ones(1,length(time));
lb = zeros(1,length(time));
Nini = 1;
U0(:,Nini) = (0.7).*rand(366,1);
Nlk = size(U0,2);
parpool('local') 
th_boh = 0:0.02:1;
Nth = length(th_boh);
Nsim = size(Lambda,2);
Uopt = zeros(Nd+1,Nsim,Nth);
Xopt = zeros(Nd+1,length(X0),Nsim,Nth);
cost_val = zeros(Nsim,Nth);
exitflags = zeros(Nsim,Nth); 
for ii = 1:Nsim
    for jj = 1:Nth
        Lam = Lambda(ii);
        param_model(4) = th_boh(jj);
        cost_functional = @(Lvect)cost_function(Lvect, param_model, param_cost, Lam, name_model);
        text = sprintf("Optimizing Lam combination %d/%d with theta = %.2f",ii,Nsim,param_model(4)); disp(text)
        options = optimoptions(@fmincon,'Algorithm','sqp', 'Display', 'iter-detailed',...
    'UseParallel', true, 'MaxIterations',150);
        [Uopt(:,ii,jj),cost_val(ii,jj),exitflags(ii,jj)] = fmincon(cost_functional,U0,[],[],[],[],lb,ub,[],options);
        Xopt(:,:,ii,jj) = ode4(name_model, time, X0, param_model, Uopt(:,ii,jj));
        clc        
    end    
end
save(title_file,'name_model','Lambda','th_boh','Uopt','Xopt','U0','cost_val','exitflags',...
    'param_model','param_cost');
delete(gcp) 