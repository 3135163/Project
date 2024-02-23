function C = cost_function(Lvect, param_model, param_cost, Lambda, name_model, idx)
Nd = param_model(18);
time = 0:Nd;
X0 = param_model(20:end);
Hbed = param_model(17);
w = param_cost(1);
VSL = param_cost(3);
Xd = param_cost(4);
p_y = param_cost(6);
cH = param_cost(7);
ts = param_cost(8);
cV = param_cost(9);
cU = param_cost(10);
X_out = ode4(name_model, time, X0, param_model, Lvect);
S_end = X_out(end,1);
H = X_out(:,6);
Q = X_out(:,7);
R1 = X_out(:,8);
dths = X_out(end,10);
N = X_out(:,11);
U=max(0,H-Hbed);
Wq = Q + min(H,Hbed) + U;
C = Lambda*( sum( w*(Lvect.*(N - Wq - ts*R1)) + Wq ) ) + ...
(1-Lambda)*( sum( cH*H + cU*U ) + (p_y*VSL + Xd)*dths );
end



