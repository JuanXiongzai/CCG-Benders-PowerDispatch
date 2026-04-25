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
% This the code for 2 nodes system test
%% coypright
% This code is part of apanj projocts in rednote
% Copyright (c) 2025 apanj
% rednote account: 9439481319
%% define varible
N_node   = 2;
node_ref = 1;
Loc_fg   = [1 1;2 1;3 2];
Loc_rg   = [1 1;2 2];
Loc_d    = [1 1;2 2];
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
P        = sdpvar(N_fg,1);
r_up     = sdpvar(N_fg,1);
r_down   = sdpvar(N_fg,1);
% power angle
sigda    = sdpvar(N_node,1); 
x        = [P;r_up;r_down;sigda];
% uncertainty variables
delta_w  = sdpvar(N_rg,1);
%%  import basic data
% transmission line
Route    = [1 ;2 ;1/0.13];      % connections of transmission lines
N_line   = size(Route,2);
Line_cap = [60];
% dispatch cost at the second stage (balancing market)
C_p      = [32;20;12];
C_ls     = [200;200];
C_ws     = zeros(N_rg,1);
C_y      = [C_p;-C_p;C_ls;C_ws;zeros(N_node,1)];
P_max    = [120;80;70];
D_max    = [110;30];
% dispatch cost at the first stage (day-ahead market)
C_rup    = [7;11;15];
C_rdown  = [5;6;14];
C_x      = [C_p;C_rup;C_rdown;zeros(N_node,1)];
% uncertainty set
w_av     = [20;25];
del_wmax = [15;20];
%% define the robust optimization problem
tem1       = Route(1,:)==[1:N_node].';
tem2       = Route(2,:)==[1:N_node].';
tem        = logical(tem1+tem2);
K          = zeros(N_node,N_node);
for i = 1:N_node
    % the conductance of the node where the power flow out
    K(i,i)                  = sum(Route(3,tem(i,:)));  
    % the conductance of the node where the power flow in
    K(i,Route(2,tem1(i,:))) = -Route(3,tem1(i,:));
    K(i,Route(1,tem2(i,:))) = -Route(3,tem2(i,:));
end
K_GL     = zeros(N_line,N_node);
for i = 1:N_line
    K_GL(i,Route(1,i)) = 1;
    K_GL(i,Route(2,i)) = -1;
    K_GL(i,:) = K_GL(i,:)*Route(3,i);
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
L             = -eye(length(y)-N_node);
L             = [L zeros(length(y)-N_node,N_node);zeros(N_line*2,length(y)-N_node) [K_GL;-K_GL]];
I             = [zeros(2*N_fg,1);-D_max;-w_av;-[Line_cap;Line_cap]];
N             = [zeros(N_fg*2+N_d,N_rg);eye(N_rg);zeros(N_line*2,N_rg)];
M             = zeros(length(y)-N_node+N_line*2,length(x));
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
g       = [-P_max;zeros(N_fg,1);-[Line_cap;Line_cap]];
% objective function at the first stage (day-ahead market)
Obj_da  = C_x'*x;
% constraints at the first stage (day-ahead market)
const_da = [F*x        == D_max-w_av
            x(end-1)   == 0
                  G*x  >= g  
           x(1:N_fg*3) >= 0
           %sum(x(N_fg+1:2*N_fg))>=26   %这行约束似乎可以删掉
                  ];
% solve the master problem
ops            = sdpsettings('solver', 'gurobi');
optimize(const_da,Obj_da,ops);
%% CCG algorithm numerical results
% initialize
LB           = -1e6;
UB           = 1e6;
epsil        = 1e-3;
M_kkt        = 1e5;
UB_iter_ccg  = [];
LB_iter_ccg  = [];
O_ccg        = [];
k_ccg        = 0;
gencons      = [];
% define variables
gencoly      = cell(10,1);
Q_w          = sdpvar(1,1);
% define master problem model
obj_MP       = Obj_da+Q_w;
% define subproblem model
delta_wup    = sdpvar(N_rg,1);
delta_wdown  = sdpvar(N_rg,1);
const_un     = [delta_wup>=0
               -delta_wup>=-del_wmax
                delta_wdown>=0
               -delta_wdown>=-del_wmax
                delta_w == delta_wup-delta_wdown
                sum((delta_wup+delta_wdown)./[15;20])<=1.4
               ];
% main function
for iter = 1:10               % assuming maximum ten iteration times
    if abs(UB-LB) <= epsil
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
           gencoly{k_ccg+1} = sdpvar(length(y),1);
           gencons          = [gencons
                               Q_w                   >= C_y'*gencoly{k_ccg+1}
                               P_cf*gencoly{k_ccg+1} == -W*delta_w_solSP-Q*x
                               L*gencoly{k_ccg+1}    >= I-M*x-N*delta_w_solSP
                               y(1:end-N_node)       >= 0
                               gencoly{k_ccg+1}(1:end-N_node) >= 0
                               gencoly{k_ccg+1}(end-1)==0
                               ];
           const_MP         = [const_da
                               gencons
                              ];
           k_ccg            = k_ccg+1;
           O_ccg            = cat(1,O_ccg,k_ccg);
        else
           % add feasibility cut
           gencoly{k_ccg+1} = sdpvar(length(y),1);
           gencons          = [gencons
                               P_cf*gencoly{k_ccg+1} == -W*delta_w_solSP-Q*x
                               L*gencoly{k_ccg+1}    >= I-M*x-N*delta_w_solSP
                               y(1:end-N_node)       >= 0
                               gencoly{k_ccg+1}(1:end-N_node) >= 0
                               gencoly{k_ccg+1}(end-1)==0
                              ];
           k_ccg            = k_ccg+1;
        end
    end
    MP_solpar    = optimize(const_MP,obj_MP,ops);
    LB           = double(obj_MP);
    LB_iter_ccg  = cat(1,LB_iter_ccg,LB);
    x_MP_sol     = double(x);
    MP_sol       = C_x'*x_MP_sol;
% subprolem
    lamb         = sdpvar(N_node,1);
    mu           = sdpvar(N_fg*2+N_d+N_rg+size(Route,2)*2,1);
    v_1          = binvar(length(mu),1);
    psi          = sdpvar(length(y)-N_node,1);
    v_2          = binvar(length(psi),1);
    ita          = sdpvar(1,1);
    const_bl_kkt = [C_y+P_cf'*lamb-L'*mu-[psi;-ita;0] == 0
                    P_cf*y                            == -W*delta_w-Q*x_MP_sol
                    y(end-1)                          == 0
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
x_robust = double(x);
%% Stochastic model
N_scen           = 200;
P_up             = sdpvar(N_fg,N_scen);
P_down           = sdpvar(N_fg,N_scen);
l_shed           = sdpvar(N_d,N_scen);
w_sp             = sdpvar(N_rg,N_scen);
Scen_weight      = 1/N_scen;
% power angle
sigbl            = sdpvar(N_node,N_scen); 
y                = [P_up;P_down;l_shed;w_sp;sigbl];
obj_stos         = Obj_da+Scen_weight*sum(C_y'*y);
delta_w_stos     = [unifrnd(-del_wmax(1),del_wmax(1),[1,N_scen]);...
                    unifrnd(-del_wmax(2),del_wmax(2),[1,N_scen])];
tem = abs(delta_w_stos(1,:))/del_wmax(1)+abs(delta_w_stos(2,:))/del_wmax(2);
while sum(tem>1.4)>0
       temreg =  tem>1.4; 
       delta_w_stos(:,temreg) = ...
                 [unifrnd(-del_wmax(1),del_wmax(1),[1,sum(temreg)]);...
                  unifrnd(-del_wmax(2),del_wmax(2),[1,sum(temreg)])];
tem = abs(delta_w_stos(1,:))/del_wmax(1)+abs(delta_w_stos(2,:))/del_wmax(2);
end
clear tem
constraints_bl_stos = [];
for i =1:N_scen
constraints_bl_stos = [constraints_bl_stos
                       P_cf*y(:,i)       == -W*delta_w_stos(:,i)-Q*x
                       L*y(:,i)          >= I-M*x-N*delta_w_stos(:,i)
                       y(1:end-N_node,i) >= 0
                       y(end-1,i)        == 0
                      ]; 
end
ops       = sdpsettings('solver', 'gurobi');
optimize([const_da constraints_bl_stos],obj_stos,ops);
x_stochas = double(x);
%% output the result of Table 1 in this paper
fprintf('Table 3\n');
fprintf(['Results for daya-ahead dispatch and reserve using the robust optimization\n' ...
    'the stochastic programmingg approach\n']);
fprintf('-----------------------------------------------------\n');
fprintf('(a) Robust  optimization  (b) Stochastic  programming\n')
fprintf('-----------------------   -------------------------\n');
fprintf(' Unit     1     2     3    Unit     1     2     3\n');
fprintf('-----------------------------------------------------\n');
for i = 1:3
    switch i
        case 1
       fprintf('P(MWh)  %3.2f  %3.2f %3.2f  P(MWh) %3.2f   %3.2f  %3.2f\n', ...
               x_robust(1), x_robust(2),x_robust(3), x_stochas(1), x_stochas(2),x_stochas(3));
        case 2
       fprintf('r+(MW)  %3.2f  %3.2f %3.2f   r+(MW) %3.2f  %3.2f   %3.2f\n', ...
               x_robust(4), x_robust(5),x_robust(6), x_stochas(4), x_stochas(5),x_stochas(6));
        case 3 
       fprintf('r-(MW)  %3.2f  %3.2f  %3.2f   r-(MW) %3.2f   %3.2f   %3.2f\n', ...
               x_robust(7), x_robust(8),x_robust(9), x_stochas(7), x_stochas(8),x_stochas(9));     
    end
end
fprintf('--------------------------------------------------------\n');