close all
clear
clc
%% C&CG algorithm for the Robust optimization
% Bo Zeng, Long Zhao,
% Solving two-stage robust optimization problems using a column-and-constraint generation method,
% Operations Research Letters,
% Volume 41, Issue 5,
% 2013,
% Pages 457-461,
% ISSN 0167-6377,
% https://doi.org/10.1016/j.orl.2013.05.003.
%% coypright
% This code is part of apanj projoct in rednote
% Copyright (c) 2025 apanj
% rednote account: 9439481319
%% C&CG algorithm numerical result
% initialize
epsil        = 1e-3;
LB           = -1e6;
UB           = 1e6;
UB_iter_ccg  = [];
LB_iter_ccg  = [];
k_ccg        = 0;
O_ccg        = [];
gencons      = [];
gencolv      = cell(10,1);
f            = [400;414;326];
a            = [18;25;20];
C            = [22 33 24;33 23 30;20 25 27];
d_par        = [206 40;274 40;220 40];
M            = 1e5;
% define variables
y            = binvar(3,1);
z            = sdpvar(3,1);
ita          = sdpvar(1);
x            = sdpvar(3,3,'full');
theta        = sdpvar(3,1);
pidual       = sdpvar(3,1);
g            = sdpvar(3,1);
v_sl         = binvar(3,1);
v_pi         = binvar(3,1);
v_x          = binvar(3,3,'full');
% define master problem
obj_MP       = f.'*y+a.'*z+ita;
cst_MPin     = [z<=800*y
                %sum(z)>=772       % this artificial calculation is removed
                z>=0
               ];
% define subproblem model
sl           = sdpvar(3,1);                  %slack variable
obj_SP       = sum(C.*x,'all')+1000*sum(sl);
d            = d_par(:,1)+d_par(:,2).*g;
cons_uncert  = [g>=0
                g<=1
                sum(g)<=1.8
                sum(g(1:2))<=1.2
               ];
% define solver options
ops          = sdpsettings('solver', 'gurobi');
% if you use other solvers, like cplex, you need to change the 'if
% condition' in the following algorithm
%% CCG algorithm main function
for iter = 1:10               % assuming maximum ten iteration times
    if UB-LB<=epsil
       break
    end
% master problem    
    if iter==1    
       const_MP = [cst_MPin
                   ita ==0
                  ];
    else
        % add optimimality cut
        % if you use other solvers, like cplex, you need to change the 'if
        % condition' in the following algorithm
        if strcmp(SP_solpar.info,'Successfully solved (GUROBI)')       
           gencolv{k_ccg+1} = sdpvar(3,3,'full');              
           gencons      = [gencons
                           ita                               >= sum(C.*gencolv{k_ccg+1},'all')
                           z-sum(gencolv{k_ccg+1},2)         >= 0
                          -u_solSP+sum(gencolv{k_ccg+1},1).' >= 0
                           reshape(gencolv{k_ccg+1},1,[])    >= 0
                          ];
           const_MP = [cst_MPin
                       gencons
                      ];
           k_ccg = k_ccg+1;
           O_ccg = cat(1,O_ccg,k_ccg);
        else
           % add feasibility cut
           gencolv{k_ccg+1} = sdpvar(3,3,'full');
           gencons          = [gencons
                               z-sum(gencolv{k_ccg+1},2)         >= 0
                              -u_solSP+sum(gencolv{k_ccg+1},1).' >= 0
                               reshape(gencolv{k_ccg+1},1,[])    >= 0
                              ];
           const_MP         = [cst_MPin
                               gencons
                               ita>=0
                              ];
           k_ccg            = k_ccg+1;
        end
    end
    MP_solpar   = optimize(const_MP,obj_MP,ops);
    LB          = double(obj_MP);
    LB_iter_ccg = cat(1,LB_iter_ccg,LB);
    MP_sol      = f.'*double(y)+a.'*double(z);
% subprolem
    %z      = [237;305;241];
    z_MPsol = double(z);
    cons_sc = [pidual           <= M*v_pi
               z_MPsol-sum(x,2) <= M*(1-v_pi)
               pidual           >= 0
               z_MPsol-sum(x,2) >= 0
               sl               <= M*v_sl
               1000+theta       <= M*(1-v_sl)
               sl               >= 0
               1000+theta       >= 0
               sum(x,1).'+sl-d  == 0
               ];
    for i = 1:3
        cons_sc = [cons_sc
                   C(i,:).'+pidual(i)+theta <= v_x(i,:).'*M
                   x(i,:).'                 <= (1-v_x(i,:).')*M
                   C(i,:).'+pidual(i)+theta >= 0
                   x(i,:).'                 >= 0
                  ];
    end
    const_SP    = [cons_sc
                   cons_uncert
                  ];        
    SP_solpar   = optimize(const_SP,-obj_SP,ops);
    if strcmp(SP_solpar.info,'Successfully solved (GUROBI)')
        UB      = min(UB,double(MP_sol)+double(obj_SP));
    end
    UB_iter_ccg = cat(1,UB_iter_ccg,UB);
    u_solSP     = double(d);
end
%% Banders-dual algorithm numerical result(comparison)
% initialize
epsil       = 1e-3;
LB          = -1e6;
UB          = 1e6;
k_ben       = 0;
O_ben       = [];
bencons     = [];
UB_iter_ben = [];
LB_iter_ben = [];
f           = [400;414;326];
a           = [18;25;20];
C           = [22 33 24;33 23 30;20 25 27];
d_par       = [206 40;274 40;220 40];
M           = 1e5;
y           = binvar(3,1);
z           = sdpvar(3,1);
ita         = sdpvar(1);
sl          = sdpvar(3,1);
theta       = sdpvar(3,1);
pidual      = sdpvar(3,1);
obj_MP      = f.'*y+a.'*z+ita;
cst_MPin    = [z      <= 800*y
               %sum(z) >= 772           % this constraint is pivotal for the feasibility of the subproblem
               z      >= 0
              ];
% uncertain set
g           = sdpvar(3,1);
d           = d_par(:,1)+d_par(:,2).*g;
cons_uncert = [g           >= 0
               g           <= 1
               sum(g)      <= 1.8
               sum(g(1:2)) <= 1.2
              ];
% solver parmeter
ops                      = sdpsettings('solver', 'gurobi','savesolveroutput',1);
ops.gurobi.InfUnbdInfo   = 1;
ops.gurobi.TuneTimeLimit = 0;    
for iter = 1:20
    if UB-LB<=epsil
       break
    end
% master problem    
    if iter==1    
       const_MP = [cst_MPin
                   ita ==0
                  ];
    else
        % add optimimality cut
        % if you use other solvers, like cplex, you need to change the 'if
        % condition' in the following algorithm
        if strcmp(SP_solpar.info,'Successfully solved (GUROBI)')      
           bencons  = [bencons
                       ita>=-z.'*double(pidual)+u_solSP.'*double(theta)
                      ];
           const_MP = [cst_MPin
                       bencons
                      ];
           k_ben    = k_ben+1;
           O_ben    = cat(1,O_ben,k_ben);
        else
           % add feasibility cut 
           exray_SP = SP_solpar.solveroutput.result.unbdray.';  
           %be careful, the feasibility cut is not tested in this problem
           bencons  = [bencons
                       exray_SP*[-z;u_solSP]<=0
                       ];
           k_ben    = k_ben+1;
        end
    end
    MP_solpar   = optimize(const_MP,obj_MP,ops);
    LB          = double(obj_MP);
    LB_iter_ben = cat(1,LB_iter_ben,LB);
    MP_sol      = f.'*double(y)+a.'*double(z);
% subprolem
    %z      = [237;305;241];
    z_MPsol = double(z);
    x       = sdpvar(3,3,'full');
    theta   = sdpvar(3,1);
    pidual  = sdpvar(3,1);
    
    v_sl = binvar(3,1);
    v_pi    = binvar(3,1);
    v_x     = binvar(3,3,'full');
    % construct subprolem model
    
    cons_sc = [pidual           <= M*v_pi
               z_MPsol-sum(x,2) <= M*(1-v_pi)
               pidual           >= 0
               z_MPsol-sum(x,2) >= 0
               sl               <= M*v_sl
               1e4-theta        <= M*(1-v_sl)
               sl               >= 0
               1e4-theta        >= 0
              -sum(x,1).'-sl+d  == 0        % I'm not sure why the equality must be written like this
               ];                           % if you write it as sum(x,1).'+sl-d==0, the algorithm cannot converge
    for i = 1:3
        cons_sc = [cons_sc
                   C(i,:).'+pidual(i)-theta <= v_x(i,:).'*M
                   x(i,:).'                 <= (1-v_x(i,:).')*M
                   C(i,:).'+pidual(i)-theta >= 0
                   x(i,:).'                 >= 0
                  ];
    end
    const_SP    = [cons_sc
                   cons_uncert
                  ];
    % 1e4 is used here as a punishment coeffience to accelerate the convergency
    obj_SP      = sum(C.*x,'all')+1e4*sum(sl);           
    SP_solpar   = optimize(const_SP,-obj_SP,ops);
    UB          = min(UB,double(MP_sol)+double(obj_SP));
    UB_iter_ben = cat(1,UB_iter_ben,UB);
    u_solSP     = double(d);
end
%% output the result of Table 1 in this paper
fprintf('Tabel 1 \nAlgorithm performance comparison\n')
fprintf('--------------------------------------------------------\n');
fprintf('  Iteration    C&CGLB     C&CGuUB     BDLB      BDUB\n');
fprintf('--------------------------------------------------------\n');
for i = 1:k_ben+1
    if k_ccg+1>=i
       fprintf('      %2d       %6.1f     %6.1f   %6.1f   %6.1f\n', ...
               i, LB_iter_ccg(i),UB_iter_ccg(i), LB_iter_ben(i),UB_iter_ben(i));
    else
       fprintf('      %2d       %6.1f     %6.1f     %6.1f   %6.1f\n', ...
               i,[] ,[] , LB_iter_ben(i),UB_iter_ben(i)); 
    end
end
%fprintf('P(A)\n');
figure
% the first three bound results are not displayed figure since its extreme
% large value makes it tricky to observe the convergency of lower and upper
% bounds
plot(3:k_ben+1,UB_iter_ben(3:end),'o-')    
hold on
plot(3:k_ben+1,LB_iter_ben(3:end),'o-')
xlabel('Iteration times')
ylabel('Bound Value')
title('Benders decomposition example')
legend('Upper bound','Lower bound')
grid on
setfig

figure
% the first three bound results are not displayed figure since its extreme
% large value makes it tricky to observe the convergency of lower and upper
% bounds
plot(2:k_ccg+1,UB_iter_ccg(2:end),'o-')    
hold on
plot(2:k_ccg+1,LB_iter_ccg(2:end),'o-')
xlabel('Iteration times')
ylabel('Bound Value')
title('CCG example')
legend('Upper bound','Lower bound')
grid on
setfig
%% conclusion
% The change in the bound of benders method is a litter bit strange here,
% and so I would like suggest you to further check KKT conditions from
% line 226 to line 235