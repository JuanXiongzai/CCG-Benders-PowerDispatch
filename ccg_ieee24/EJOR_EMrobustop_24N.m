clear
close all
clc
%% literature
% Marco Zugno, Antonio J. Conejo,
% A robust optimization approach to energy and reserve dispatch in electricity markets,
% European Journal of Operational Research,
% Volume 247, Issue 2,
% 2015,
% Pages 659-671,
% https://doi.org/10.1016/j.ejor.2015.05.081.
% This the code for 24 nodes system test
%% coypright
% This code is part of apanj projocts in rednote
% Copyright (c) 2025 apanj
% rednote account: 9439481319
%% import data
TL        = readtable('IEEE24RTSDataROP.xlsx','Sheet','TL','VariableNamingRule','preserve');
FFgen     = readtable('IEEE24RTSDataROP.xlsx','Sheet','FFgen','VariableNamingRule','preserve');
RESgen    = readtable('IEEE24RTSDataROP.xlsx','Sheet','RESgen','VariableNamingRule','preserve');
LoadDis   = readtable('IEEE24RTSDataROP.xlsx','Sheet','LoadDis','VariableNamingRule','preserve');
FFgenCost = readtable('IEEE24RTSDataROP.xlsx','Sheet','FFgenCost','VariableNamingRule','preserve');
%% define varible
N_node   = 24;
node_ref = 1;
Loc_fg   = [FFgen.Unit FFgen.Node];
Loc_rg   = [RESgen.Unit RESgen.Node];
Loc_d    = [LoadDis.Load LoadDis.Node];
N_fg     = length(Loc_fg);     % fossil-fuel generator
N_rg     = length(Loc_rg);     % renewable generator
N_d      = length(Loc_d);
% decision variables at the second stage (balancing market)
P_up     = sdpvar(N_fg,1);
P_down   = sdpvar(N_fg,1);
l_shed   = sdpvar(N_d,1);
w_sp     = sdpvar(N_rg,1);
% power angle
sigbl    = sdpvar(N_node,1); 
y        = [P_up;P_down;l_shed;w_sp;sigbl];
% decision variables at the first stage (day-ahead market)
P        = sdpvar(N_fg,1,'full');
r_up     = sdpvar(N_fg,1,'full');
r_down   = sdpvar(N_fg,1,'full');
% power angle
sigda    = sdpvar(N_node,1,'full'); 
x        = [P;r_up;r_down;sigda];
% uncertainty variables
delta_w  = sdpvar(N_rg,1);
%%  import basic data
% transmission line
Route             = [TL.From TL.To 1./TL.Reactance].';      % connections of transmission lines
N_line            = size(Route,2);
Line_cap          = TL.Capacity;
% dispatch cost at the second stage (balancing market)
C_regup           = FFgenCost.Cregup;
C_regdown         = FFgenCost.Cregdown;
C_ls              = 200*ones(N_d,1);
C_ws              = zeros(N_rg,1);
C_y               = [C_regup;-C_regdown;C_ls;C_ws;zeros(N_node,1)];
P_max             = FFgen.Pmax;
D_max             = zeros(N_node,1);
D_max(LoadDis.Node) = LoadDis.percentage*2850;
% dispatch cost at the first stage (day-ahead market)
C_p               = FFgenCost.C;
C_rup             = FFgenCost.Cresup;
C_rdown           = FFgenCost.Cresdown;
C_x               = [C_p;C_rup;C_rdown;zeros(N_node,1)];
% uncertainty set
load avgwpower.mat  % average windpower output based on 2023 wind speed data
w_av              = zeros(N_node,1);
w_av(Loc_rg(:,2)) = avgwpower;
wratepower        = [420;315;455;245;175;420];
delta_wmax        = max(wratepower-avgwpower,avgwpower);
%% define the robust optimization problem
tem1       = Route(1,:)== [1:N_node].';
tem2       = Route(2,:)== [1:N_node].';
tem        = logical(tem1+tem2);
K          = zeros(N_node,N_node);
for i = 1:N_node
    % the conductance of the node where the power flow out
    K(i,i)                  = sum(Route(3,tem(i,:)));  
    % the conductance of the node where the power flow in
    K(i,Route(2,tem1(i,:))) = -Route(3,tem1(i,:));
    K(i,Route(1,tem2(i,:))) = -Route(3,tem2(i,:));
end
K_GL  = zeros(N_line,N_node);
for i = 1:N_line
    K_GL(i,Route(1,i)) = 1;
    K_GL(i,Route(2,i)) = -1;
    K_GL(i,:)          = K_GL(i,:)*Route(3,i);
end
% objective function at the second stage (balancing market)
Obj_bl     = C_y'*y;
% coefficent of equality constraints
P_cf_pup   = (Loc_fg(:,2)==[1:N_node])';          % cf refers to coefficent 
P_cf_pdown = -P_cf_pup;
P_cf_wp    = -(Loc_rg(:,2)==[1:N_node])';
P_cf_ls    = (Loc_d(:,2)==[1:N_node])';
P_cf       = [P_cf_pup P_cf_pdown P_cf_ls P_cf_wp -K]; 
W          = -P_cf_wp;                            % fluctuation in wind power output
Q          = [zeros(N_node,N_fg*3) K];
%K(:,node_ref) = 0;
% coefficent of inequality constraints
L          = -eye(length(y)-N_node);
L          = [L zeros(length(y)-N_node,N_node);zeros(N_line*2,length(y)-N_node) [K_GL;-K_GL]];
I          = [zeros(2*N_fg,1);-LoadDis.percentage*2850;-avgwpower;-[Line_cap;Line_cap]];
N          = [zeros(N_fg*2+N_d,N_rg);eye(N_rg);zeros(N_line*2,N_rg)];
M          = zeros(length(y)-N_node+N_line*2,length(x));
M(1:N_fg*2,N_fg+1:N_fg*3) = eye(N_fg*2);
clear tem1 tem2 tem
%% define the master problem
% coefficent of equality constraints
F       = [P_cf_pup zeros(N_node,N_fg*2) -K];
G_p     = [-eye(N_fg);eye(N_fg)];
G_rup   = [-eye(N_fg);zeros(N_fg,N_fg)];
G_rdown = [zeros(N_fg,N_fg);-eye(N_fg);];
% coefficent of inequality constraints
G       = [[G_p G_rup G_rdown;zeros(N_line*2,N_fg*3)] [zeros(N_fg*2,N_node);K_GL;-K_GL]];
G       = [G;[zeros(N_fg*2,N_fg)  -eye(N_fg*2)  zeros(N_fg*2,N_node)]];
g       = [-P_max;zeros(N_fg,1);-[Line_cap;Line_cap]];
g       = [g;-FFgen.Resup;-FFgen.Resdown];                   %reserve capacity
% objective function at the first stage (day-ahead market)
Obj_da  = C_x'*x;
% constraints at the first stage (day-ahead market)
const_da = [F*x         == D_max-w_av
            x(N_fg*3+1) == 0                % 1 node is the reference node
                  G*x   >= g  
            x(1:N_fg*3) >= 0           
           ]; 
% solve the master problem
ops            = sdpsettings('solver', 'gurobi');
%saveampl(const_da, Obj_da,'mymodelDA')
optimize(const_da,Obj_da,ops);
%xtest = double(x);
% test balancing market
%sig_MP_sol  = double(sigda);
%x_MP_sol    = double(x);
%delta_w_test = [-162.94;-3.959;-172.650;-26.682;-40.808;-228.701];
%constraints_bl_test = [P_cf*y               == -W*delta_w_test-Q*x_MP_sol
%                       L*y                  >= I-M*x_MP_sol-N*delta_w_test
%                       y(1:end-N_node)      >= 0
%                       y(2*N_fg+N_d+N_rg+1) == 0
%                       ]; 
%constraints_SP = [constraints_bl_test
%                  const_un
%                 ];
%ops            = sdpsettings('solver', 'gurobi');
%SP_solpar      = optimize(constraints_SP,Obj_bl,ops);
%ytest = double(y);
%P_cf*ytest-(-W*delta_w_test-Q*x_MP_sol)
%% CCG algorithm numerical results
% initialize
LB           = -1e10;
UB           = 1e10;
epsil        = 1e-3;
M_kkt        = 1e5;
UB_iter_ccg  = [];
LB_iter_ccg  = [];
O_ccg        = [];
k_ccg        = 0;
gencons      = [];
% define variables
gencoly      = cell(10,1);                  %assume the maximum iterations
Q_w          = sdpvar(1,1);                 %optimization function in the second stage
% define master problem model
obj_MP       = Obj_da+Q_w;
% define subproblem model
delta_wup    = sdpvar(N_rg,1);
delta_wdown  = sdpvar(N_rg,1);
const_un     = [delta_wup>=0
                delta_wup<=delta_wmax
                delta_wdown>=0
               delta_wdown<=delta_wmax
               delta_w == delta_wup-delta_wdown
               avgwpower+delta_w>=0
               sum((delta_wup+delta_wdown)./delta_wmax)<=4.71
               ];

% main function
for iter = 1:10               % assuming maximum ten iteration times
    if abs(UB-LB)<=epsil
       break
    end
    % master problem    
    if iter==1    
       const_MP = [const_da
                   Q_w==0
                  ];
    else
        % add optimimality cut
        if strcmp(SP_solpar.info,'Successfully solved (GUROBI)')       
           gencoly{k_ccg+1}   = sdpvar(length(y),1);
           gencons      = [gencons
                           Q_w                                 >= C_y'*gencoly{k_ccg+1}
                           P_cf*gencoly{k_ccg+1}               == -W*delta_w_solSP-Q*x
                           L*gencoly{k_ccg+1}                  >= I-M*x-N*delta_w_solSP
                           y(1:end-N_node)                     >= 0
                           gencoly{k_ccg+1}(1:end-N_node)      >= 0
                           gencoly{k_ccg+1}(2*N_fg+N_d+N_rg+1) == 0
                           ];
           const_MP = [const_da
                       gencons
                      ];
           k_ccg = k_ccg+1;
           O_ccg = cat(1,O_ccg,k_ccg);
        else
           % add feasibility cut
           gencoly{k_ccg+1}   = sdpvar(length(y),1);
           gencons            = [gencons
                                 P_cf*gencoly{k_ccg+1}               == -W*delta_w_solSP-Q*x
                                 L*gencoly{k_ccg+1}                  >= I-M*x-N*delta_w_solSP
                                 y(1:end-N_node)                     >= 0
                                 gencoly{k_ccg+1}(1:end-N_node)      >= 0
                                 gencoly{k_ccg+1}(2*N_fg+N_d+N_rg+1) == 0
                               ];
           k_ccg = k_ccg+1;
        end
    end
    MP_solpar   = optimize(const_MP,obj_MP,ops);
    LB          = double(obj_MP);
    LB_iter_ccg = cat(1,LB_iter_ccg,LB);
    x_MP_sol    = double(x);
    MP_sol      = C_x'*x_MP_sol;
% subprolem
    lamb         = sdpvar(N_node,1);
    mu           = sdpvar(N_fg*2+N_d+N_rg+N_line*2,1);
    v_1          = binvar(length(mu),1);
    psi          = sdpvar(length(y)-N_node,1);
    v_2          = binvar(length(psi),1);
    ita          = sdpvar(1,1);
    const_bl_kkt = [C_y+P_cf'*lamb-L'*mu-[psi;-ita;zeros(N_node-1,1)] == 0
                    P_cf*y                            == -W*delta_w-Q*x_MP_sol
                    y(2*N_fg+N_d+N_rg+1)                   == 0
                    L*y-(I-M*x_MP_sol-N*delta_w)      <= M_kkt*(1-v_1)
                    L*y-(I-M*x_MP_sol-N*delta_w)      >= 0
                    mu                                <= M_kkt*v_1
                    mu                                >= 0
                    y(1:end-N_node)                   <= M_kkt*(1-v_2)
                    y(1:end-N_node)                   >= 0
                    psi                               <= M_kkt*v_2
                    psi                               >= 0 
                    ];
    const_SP     = [const_bl_kkt
                    const_un
                    ];
   SP_solpar     = optimize(const_SP,-Obj_bl,ops);
   delta_w_solSP = double(delta_w);
   UB            = min(UB,MP_sol+double(Obj_bl));
   UB_iter_ccg   = cat(1,UB_iter_ccg,UB);
end                                                                                   
x_robust         = double(x);
y_robust         = double(y);
%% Stochastic model
% The windpower generated from 2024 wind speed data is used to solve the
% stochastic model
N_scen         = 500;                          % number of scenario
P_up           = sdpvar(N_fg,N_scen,'full');
P_down         = sdpvar(N_fg,N_scen,'full');
l_shed         = sdpvar(N_d,N_scen,'full');
w_sp           = sdpvar(N_rg,N_scen,'full');
Scen_weight    = 1/N_scen;
% power angle
sigbl          = sdpvar(N_node,N_scen,'full'); 
y_stocs        = [P_up;P_down;l_shed;w_sp;sigbl];

obj_stos       = Obj_da+Scen_weight*sum(C_y'*y_stocs);
delta_w_stos   = struct2array(load('wpower_unit2023.mat'))-avgwpower;
%delta_w_stos   = zeros(N_rg,N_scen);
%for i = 1:N_rg
%    delta_w_stos(i,:) = unifrnd(0,wratepower(i),[1,N_scen])-avgwpower(i);
%end
%tem = abs(sum(delta_w_stos./delta_wmax));
%while sum(tem>2.71)>0
%       temreg =  tem>2.71; 
%       for i = 1:N_rg
%           delta_w_stos(i,temreg) = unifrnd(0,wratepower(i),[1,sum(temreg)])-avgwpower(i);
%       end
%tem = abs(sum(delta_w_stos./delta_wmax));
%end
%clear tem
constraints_bl_stos = [];
for i =1:N_scen
constraints_bl_stos = [constraints_bl_stos
                       P_cf*y_stocs(:,i)            == -W*delta_w_stos(:,i)-Q*x
                       L*y_stocs(:,i)               >= I-M*x-N*delta_w_stos(:,i)
                       y_stocs(1:end-N_node,i)      >= 0
                       y_stocs(2*N_fg+N_d+N_rg+1,i) == 0
                       ]; 
end
optimize([const_da constraints_bl_stos],obj_stos,ops);
x_stochas       = double(x);
%[model,recoverymodel] = export([const_da constraints_bl_stos], obj_stos, sdpsettings('solver','GUROBI'));
%iis                   = gurobi_iis(model);
%const_test            = [const_da constraints_bl_stos];
%find(iis.Arows==1)
%% Test in the out-of-sample
alpha_out      = 0.95;
delta_w_test   = struct2array(load('wpower_unit2024.mat'))-avgwpower;
P_up           = sdpvar(N_fg,1,'full');
P_down         = sdpvar(N_fg,1,'full');
l_shed         = sdpvar(N_d,1,'full');
w_sp           = sdpvar(N_rg,1,'full');
% power angle
sigbl          = sdpvar(N_node,1,'full'); 
y_test         = [P_up;P_down;l_shed;w_sp;sigbl];
Obj_bltest     = C_y'*y_test;
ops            = sdpsettings('solver', 'gurobi');
% Stochastic model
Balance_cost   = zeros(N_scen,1);
Energy_redisp  = zeros(N_scen,1);
Load_shed      = zeros(N_scen,1);
for i = 1:N_scen
    constraints_bltest = [P_cf*y_test               == -W*delta_w_test(:,i)-Q*x_stochas
                          L*y_test                  >= I-M*x_stochas-N*delta_w_test(:,i)
                          y_test(1:end-N_node)      >= 0
                          y_test(2*N_fg+N_d+N_rg+1) == 0
                         ]; 
    % solve the master problem    
    ofs_sol           = optimize(constraints_bltest,Obj_bltest,ops);
    y_test_sol        = double(y_test);
    Balance_cost(i)   = double(Obj_bltest);
    Energy_redisp(i)  = C_y(1:2*N_fg).'*y_test_sol(1:2*N_fg);
    Load_shed(i)      = C_y(2*N_fg+1:2*N_fg+N_d).'*y_test_sol(2*N_fg+1:2*N_fg+N_d);     
end
% Calculate CVaR results and VaR results (the worst scenario)
CVaR_Baldata_stos = maxk(Balance_cost,ceil((1-alpha_out)*500));
CVaR_Bal95_stos   = sum(CVaR_Baldata_stos)/(500*(1-alpha_out));
VaR_Bal100_stos   = max(CVaR_Baldata_stos);
CVaR_Engdata_stos = maxk(Energy_redisp,ceil((1-alpha_out)*500));
CVaR_Eng95_stos   = sum(CVaR_Engdata_stos)/(500*(1-alpha_out));
VaR_Eng100_stos   = max(CVaR_Engdata_stos);
CVaR_Lsheddata_stos = maxk(Load_shed,ceil((1-alpha_out)*500));
CVaR_Lshed95_stos = sum(CVaR_Lsheddata_stos)/(500*(1-alpha_out));
VaR_Lshed100_stos = max(CVaR_Lsheddata_stos);
clear CVaR_Baldata_stos CVaR_Engdata_stos CVaR_Lsheddata_stos

Avg_Bal_stos      = mean(Balance_cost);
Avg_Eng_stos      = mean(Energy_redisp);
Avg_Lshed_stos    = mean(Load_shed);
Engdisp_stos      = C_x(1:N_fg).'*x_stochas(1:N_fg);
Upreserve_stos    = C_x(N_fg+1:2*N_fg).'*x_stochas(N_fg+1:2*N_fg);  
Downreserve_stos  = C_x(2*N_fg+1:3*N_fg).'*x_stochas(2*N_fg+1:3*N_fg);  
DACost_stos       = C_x.'*x_stochas;

Avg_Tot_stos      = DACost_stos+Avg_Bal_stos;
CVaR_Tot95_stos   = DACost_stos+CVaR_Bal95_stos;
VaR_Tot100_stos   = DACost_stos+VaR_Bal100_stos;
% Robust model
Balance_cost      = zeros(N_scen,1);
Energy_redisp     = zeros(N_scen,1);
Load_shed         = zeros(N_scen,1);
for i = 1:N_scen
    constraints_bltest = [P_cf*y_test               == -W*delta_w_test(:,i)-Q*x_robust
                          L*y_test                  >= I-M*x_robust-N*delta_w_test(:,i)
                          y_test(1:end-N_node)      >= 0
                          y_test(2*N_fg+N_d+N_rg+1) == 0
                         ]; 
    % solve the master problem    
    ofs_sol           = optimize(constraints_bltest,Obj_bltest,ops);
    y_test_sol        = double(y_test);
    Balance_cost(i)   = double(Obj_bltest);
    Energy_redisp(i)  = C_y(1:2*N_fg).'*y_test_sol(1:2*N_fg);
    Load_shed(i)      = C_y(2*N_fg+1:2*N_fg+N_d).'*y_test_sol(2*N_fg+1:2*N_fg+N_d);     
end
% Calculate CVaR results and VaR results (the worst scenario)
CVaR_Baldata_rob = maxk(Balance_cost,ceil((1-alpha_out)*500));
CVaR_Bal95_rob   = sum(CVaR_Baldata_rob)/(500*(1-alpha_out));
VaR_Bal100_rob   = max(CVaR_Baldata_rob);
CVaR_Engdata_rob = maxk(Energy_redisp,ceil((1-alpha_out)*500));
CVaR_Eng95_rob   = sum(CVaR_Engdata_rob)/(500*(1-alpha_out));
VaR_Eng100_rob   = max(CVaR_Engdata_rob);
CVaR_Lsheddata_rob = maxk(Load_shed,ceil((1-alpha_out)*500));
CVaR_Lshed95_rob = sum(CVaR_Lsheddata_rob)/(500*(1-alpha_out));
VaR_Lshed100_rob = max(CVaR_Lsheddata_rob);
clear CVaR_Baldata_rob CVaR_Engdata_rob CVaR_Lsheddata_rob

Avg_Bal_rob      = mean(Balance_cost);
Avg_Eng_rob      = mean(Energy_redisp);
Avg_Lshed_rob    = mean(Load_shed);
Engdisp_rob      = C_x(1:N_fg).'*x_robust(1:N_fg);
Upreserve_rob    = C_x(N_fg+1:2*N_fg).'*x_robust(N_fg+1:2*N_fg);  
Downreserve_rob  = C_x(2*N_fg+1:3*N_fg).'*x_robust(2*N_fg+1:3*N_fg);
DACost_rob       = C_x.'*x_robust;

Avg_Tot_rob      = DACost_rob+Avg_Bal_rob;
CVaR_Tot95_rob   = DACost_rob+CVaR_Bal95_rob;
VaR_Tot100_rob   = DACost_rob+VaR_Bal100_rob;
%% output the result of Table 1 in this paper
fprintf('Table 6\n');
fprintf('Comparison of system cost: robust optimization, stochastic programming\n');
fprintf('----------------------------------------------------------------------------------------\n');
fprintf('Cost(dollar)             (a) Robust  optimization       (b) Stochastic  programming\n')
fprintf('----------------------------------------------------------------------------------------\n');
fprintf(' Energydispatch              %8.2f                       %8.2f \n',Engdisp_rob,Engdisp_stos);
fprintf(' Upward reserve              %8.2f                       %8.2f \n',Upreserve_rob,Upreserve_stos);
fprintf(' downward reserve            %8.2f                       %8.2f \n',Downreserve_rob,Downreserve_stos);
fprintf(' Day-ahead market            %8.2f                       %8.2f \n',DACost_rob,DACost_stos);
fprintf('\n');
fprintf('                           Mean     CVaR95    VaR100     Mean     CVaR95     VaR100  \n');
fprintf(' Energy redispatch         %7.2f  %7.2f  %7.2f   %7.2f   %7.2f    %7.2f\n',Avg_Eng_rob,CVaR_Eng95_rob,VaR_Eng100_rob,...
        Avg_Eng_stos,CVaR_Eng95_stos,VaR_Eng100_stos);
fprintf(' Load-shedding             %7.2f  %7.2f  %7.2f   %7.2f  %7.2f   %7.2f\n',Avg_Lshed_rob,CVaR_Lshed95_rob,VaR_Lshed100_rob,...
        Avg_Lshed_stos,CVaR_Lshed95_stos,VaR_Lshed100_stos);
fprintf(' Balancing market          %7.2f  %7.2f  %7.2f   %7.2f  %7.2f   %7.2f\n',Avg_Bal_rob,CVaR_Bal95_rob,VaR_Bal100_rob,...
        Avg_Bal_stos,CVaR_Bal95_stos,VaR_Bal100_stos);
fprintf(' Electricty market        %7.2f  %7.2f  %7.2f  %7.2f  %7.2f   %7.2f\n',Avg_Tot_rob,CVaR_Tot95_rob,VaR_Tot100_rob,...
        Avg_Tot_stos,CVaR_Tot95_stos,VaR_Tot100_stos);
fprintf('----------------------------------------------------------------------------------------\n');
%% CVaR sensitive analysis
[x_upres_stos,x_downres_stos] = CVaRsensitive(N_scen,N_fg,Obj_da,x,C_y,y_stocs,const_da,constraints_bl_stos);
% I think it is unnecessary to 'clear' the down reserve in this problem
% since there is no punishment in the balancing market, i.e, the wind spillrage
% is free of cost.
x_upres_rob = sum(x_robust(N_fg+1:2*N_fg));
alpha        = 0:0.1:0.9;
alpha        = [alpha 0.99];
figure
plot(alpha,x_upres_stos)
hold on
grid on
line([0 1],[x_upres_rob x_upres_rob],'Color','red','LineStyle','--')
ylim([0 600])
xlabel('$\alpha$','Interpreter','latex')
ylabel('Up reserve capacity(MW)')
legend('Stochastic model','Robust model')
title({'The CVaR sensitive analysis for','up-reserve capacity cleared in DA market'})