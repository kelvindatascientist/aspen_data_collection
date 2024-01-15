function [CBM_cond,CBM0_cond,CBM_reb,CBM0_reb]=Size_RadFRac_condreb(simulation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Condensador
Type_A = 'Heat Exchanger';
Type_B = 'Heat Exchanger Floating Head';
% o condensador é um trocador de casco e tubo feito de aço carbono
% refrigeração água a 4bar e 30ºC, saindo do trocador a 40ºC
U_cond = 765;  %W/m2.K, tabela do Campbell (2011) Figure 13.2
T_supply = 20; %oC, considerando 4 bar
T_return = 30; %oC , pode ser 50oC precisa de ref.
Q_cond = simulation{1};  % W , Obtido da simulação
T_top = simulation{2} % Application.Tree.FindNode("\Data\Blocks\DC\Output\B_TEMP\1")
P_cond= simulation{3}; %Application.Tree.FindNode("\Data\Blocks\DC\Input\PRES1")
LMTD(T_supply,T_return,T_top,T_top)
A_req_cond = -Q_cond/(U_cond*LMTD(T_supply,T_return,T_top,T_top));

A_cond_min  = 10;   
A_cond_max  = 1000;

%%% Calculating the minimum number of units
N_HeatEx_cond = ceil(A_req_cond/A_cond_max);    	

%%% Sizing
A_Prelim = A_req_cond/N_HeatEx_cond;              
A_cond = max(A_Prelim,A_cond_min)
% Sizing = {Cap_Unit, P, D_Col,N_tray, N_Spares,Type_Material,Type_Material_Tube,Type_A,Type_B}
Sizing={A_cond, P_cond,[],[],[],'CS','CS',Type_A,Type_B}
[CBM_cond,CBM0_cond]=CapitalCost(Sizing);
CBM_cond=CBM_cond*N_HeatEx_cond
CBM0_cond=CBM0_cond*N_HeatEx_cond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%reboiler
Type_A = 'Heat Exchanger';
Type_B = 'Heat Exchanger Kettle Reboiler';
% refervedor é um kettle reboiler de aço carbono
% aquecimento no refeverdor, utiliza-se, quando possível, vapor saturado de 
% baixa pressão - 6 bar (160ºC), ou de média pressão - 10 bar (184ºC)
U_cond = 910;  %W/m2.K, tabela do Campbell (2011) Figure 13.2
T_supply = 160; %oC, considerando vapor fraction 1
T_return = 160; %oC , considerando vapor fraction 0
Q_reb = simulation{4}  % W , Obtido da simulação
T_bot_L = simulation{5} % Application.Tree.FindNode("\Data\Blocks\DC\Output\B_TEMP\1")
T_bot_V = simulation{6}
P_bot=simulation{7}
LMTD(T_bot_L,T_bot_V,T_supply,T_return)
A_req_reb = Q_reb/(U_cond*LMTD(T_bot_L,T_bot_L,T_supply,T_return))

A_reb_min  = 10;   
A_reb_max  = 1000;

%%% Calculating the minimum number of units
N_HeatEx_reb = ceil(A_req_reb/A_reb_max)   	

%%% Sizing
A_Prelim = A_req_reb/N_HeatEx_reb            
A_reb = max(A_Prelim,A_reb_min)
% Sizing = {Cap_Unit, P, D_Col,N_tray, N_Spares,Type_Material,Type_Material_Tube,Type_A,Type_B}
Sizing={A_reb, P_bot,[],[],[],'CS','CS',Type_A,Type_B}
[CBM_reb,CBM0_reb]=CapitalCost(Sizing);
CBM_reb=CBM_reb*N_HeatEx_reb
CBM0_reb=CBM0_reb*N_HeatEx_reb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function log_mean = LMTD(T_CI,T_CO,T_HI,T_HO)
if (T_HI - T_CO) - (T_HO - T_CI) == 0;
    log_mean = abs(T_HI - T_CO);
else
    log_mean = ((T_HI - T_CO) - (T_HO - T_CI))/log((T_HI - T_CO)/(T_HO - T_CI));
end
end

end