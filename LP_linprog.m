%%%%%%%%%%%%%%%%%%%%%%%% Created by Mohsen Lotfi %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2022/09/15 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Linear programming problem (LPP): ==> linprog
format short
clear all
clc
%% read the data
filename = 'data.xlsx';
sheet = 2;
xlRange_1 = 'A2:A97';
xlRange_2 = 'B2:B97';
demand = xlsread(filename,sheet,xlRange_1);
price_el = xlsread(filename,sheet,xlRange_2);

%%
% Constant
%x0 = [0.5*ones(1,192),zeros(1,95)];

P_gas  = 30;

eta_gas  = 0.85;
eta_el = 0.95;

% C[1] = C[97] =0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% g[i] E[i] C[i] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective
F = [30*ones(1,96),price_el',zeros(1,95)];
    
% inequality constraints
A_ineq = [];
for i=0:93
    A_aux_ineq = [zeros(1,192+i),1,-1,zeros(1,93-i);zeros(1,192+i),-1,1,zeros(1,93-i)];
    A_ineq     = [A_ineq;A_aux_ineq];
end
b_ineq = 2*ones(1,188);

% equality constraints
A_eq = [];
for i=1:94
    A_aux_eq = [zeros(1,i),eta_gas,zeros(1,95),eta_el,zeros(1,94),1,-1,zeros(1,94-i)];
    A_eq     = [A_eq;A_aux_eq];
end
row_1  = [eta_gas,zeros(1,95),eta_el,zeros(1,95),-1,zeros(1,94)];
row_96 = [zeros(1,95),eta_gas,zeros(1,95),eta_el,zeros(1,94),1];
A_eq = [row_1;A_eq;row_96];
b_eq = demand;

% bounds
lb = zeros(1,287);
ub = [3*ones(1,96),2*ones(1,96),6*ones(1,95)]; %C[2]=C[96]=2

options = optimoptions('linprog','Algorithm','dual-simplex');

% linprog function
X = linprog(F,A_ineq,b_ineq,A_eq,b_eq,lb,ub,options);

% optimum price
optimum_price = 30*ones(1,96)*X(1:96)+price_el'*X(97:192)

%% Plot
figure(1)
plot(1:96,X(1:96),'-','MarkerSize',10)
hold on
plot(1:96,X(97:192),'-','MarkerSize',10)
legend('\it Gas boiler production','\it Electrode boiler production')
xlabel('Time [1h]','Interpreter','LaTex')
ylabel('Thermal energy','Interpreter','LaTex')
title('Thermal energy without optimization','Interpreter','LaTex')
set(gca,'FontName','Times New Roman','FontSize',12,'fontWeight','bold');
grid on

figure(2)
plot(1:96,eta_gas*X(1:96)+eta_el*X(97:192),1:96,demand,'-','MarkerSize',10)
hold on
plot(1:96,demand-eta_gas*X(1:96)-eta_el*X(97:192),'-','MarkerSize',10)
legend('\it Mixed Production','\it Thermal energy demand','\it Thermal energy storage')
xlabel('Time','Interpreter','LaTex')
ylabel('Thermal energy','Interpreter','LaTex')
title('Thermal energy with optimization','Interpreter','LaTex')
set(gca,'FontName','Times New Roman','FontSize',12,'fontWeight','bold');
grid on

