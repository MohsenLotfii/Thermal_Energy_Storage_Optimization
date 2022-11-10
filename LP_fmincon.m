%%%%%%%%%%%%%%%%%%%%%%%% Created by Mohsen Lotfi %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2022/09/12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Linear programming problem (LPP)==> fmincon
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

%% Solve fmincon
% Initial point in the middle of the region
x0 = [0.75*ones(1,192),zeros(1,95)];

% Constants:
P_gas   = 30;
eta_gas = 0.85;
eta_el  = 0.95;
% C[1]=C[97]=0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% g[i] E[i] C[i] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function: is a simple algebraic function of two variables
F = @(x)P_gas*ones(1,96)*x(1:96)'+price_el'*x(97:192)';

% inequality constraints
A_ineq = [];
for j=0:93
    A_aux_ineq = [zeros(1,192+j),1,-1,zeros(1,93-j);zeros(1,192+j),-1,1,zeros(1,93-j)];
    A_ineq     = [A_ineq;A_aux_ineq];
end
b_ineq = 2*ones(1,188); % for C[i] ==> 94*2=188

% equality constraints
A_eq = [];
row_1  = [eta_gas,zeros(1,95),eta_el,zeros(1,95),-1,zeros(1,94)];
for i=1:94
    A_aux_eq = [zeros(1,i),eta_gas,zeros(1,95),eta_el,zeros(1,94),1,-1,zeros(1,94-i)];
    A_eq     = [A_eq;A_aux_eq];
end
row_96 = [zeros(1,95),eta_gas,zeros(1,95),eta_el,zeros(1,94),1];
A_eq = [row_1;A_eq;row_96];
b_eq = demand;

% bounds
lb = zeros(1,287); % 96(for g[i])+96(for E[i])+95(for C[i])=287
ub = [3*ones(1,96),2*ones(1,96),6*ones(1,95)];

% fmincon function: x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub)
X = fmincon(F,x0,A_ineq,b_ineq,A_eq,b_eq,lb,ub);

% optimum price
optimum_price = P_gas*ones(1,96)*X(1:96)'+price_el'*X(97:192)'

%% Plot
figure(1)
plot(1:96,X(1:96),'-','MarkerSize',10)
hold on
plot(1:96,X(97:192),'-','MarkerSize',10)
legend('\it Gas boiler production','\it Electrode boiler production')
xlabel('Time [1h]','Interpreter','LaTex')
ylabel('Energy production [MW]','Interpreter','LaTex')
title('Thermal energy with optimization','Interpreter','LaTex')
set(gca,'FontName','Times New Roman','FontSize',12,'fontWeight','bold');
grid on

figure(2)
plot(1:96,eta_gas*X(1:96)+eta_el*X(97:192),1:96,demand,'-','MarkerSize',10)
hold on
plot(1:96,demand'-eta_gas*X(1:96)-eta_el*X(97:192),'-','MarkerSize',10)
legend('\it Mixed Production','\it Thermal energy demand','\it Thermal energy storage')
xlabel('Time','Interpreter','LaTex')
ylabel('Energy production [MW]','Interpreter','LaTex')
title('Thermal energy with optimization','Interpreter','LaTex')
set(gca,'FontName','Times New Roman','FontSize',12,'fontWeight','bold');
grid on

