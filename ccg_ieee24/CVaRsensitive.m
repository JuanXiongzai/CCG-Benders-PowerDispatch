function [x_upres,x_downres] = CVaRsensitive(N_scen,N_fg,Obj_da,x,C_y,y_stocs,const_da,constraints_bl_stos)
Scen_weight      = 1/N_scen;
obj_stos         = Obj_da+Scen_weight*sum(C_y'*y_stocs);
beta             = 50;
zeta             = sdpvar(1,1);
lamb             = sdpvar(N_scen,1);
ops              = sdpsettings('solver', 'gurobi');
x_upres          = zeros(10,1);
x_downres        = zeros(10,1);
constraints_CvaR = [lamb>=0
                    (Obj_da+(C_y'*y_stocs).')-zeta-lamb<=0;
                   ];
for i = 1:11
    if i<11
       alpha     = 0.1*(i-1);
    else
       alpha     = 0.1*(i-1)-0.01;
    end
    CvaR         = Scen_weight*sum(lamb)/(1-alpha)+zeta;
    sol          = optimize([const_da constraints_bl_stos constraints_CvaR],...
                             obj_stos+beta*CvaR,ops);
    x_stochas    = double(x);
    x_upres(i)   = sum(x_stochas(N_fg+1:2*N_fg));
    x_downres(i) = sum(x_stochas(2*N_fg+1:3*N_fg));
end