function [Sizing,N_HeatEx_cond]=Size_coller(simulation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:
% simulation(1) = Q
% simulation(2) = 
%Condensador
Type_A = 'Heat Exchanger';
Type_B = 'Heat Exchanger Floating Head';
utility_cond = simulation{5};
if (strcmp(utility_cond,'WATER'))
    % o condensador é um trocador de casco e tubo feito de aço carbono
    % refrigeração água a 4bar e 30ºC, saindo do trocador a 40ºC
    U_cond = 765;  %W/m2.K, tabela do Campbell (2011) Figure 13.2
    T_supply = 20; %oC, considerando 4 bar
    T_return = 30; %oC , pode ser 50oC precisa de ref.
elseif (strcmp(utility_cond,'FREON'))
        U_cond = 765;  %W/m2.K,
        T_supply = -25; %oC,
        T_return = -25; %oC ,
end
Q_cond = simulation{1};  % W , Obtido da simulação
T_in = simulation{2};
T_out = simulation{3};
P_cond= simulation{4};
DLMTD = LMTD(T_supply,T_return,T_in,T_out)
A_req_cond = -Q_cond/(U_cond*LMTD(T_supply,T_return,T_in,T_out));

A_cond_min  = 10;   
A_cond_max  = 1000;

%%% Calculating the minimum number of units
N_HeatEx_cond = ceil(A_req_cond/A_cond_max)    	

%%% Sizing
A_Prelim = A_req_cond/N_HeatEx_cond;              
A_cond = max(A_Prelim,A_cond_min)

% Sizing = {Cap_Unit, P, D_Col,N_tray, N_Spares,Type_Material,Type_Material_Tube,Type_A,Type_B}
Type_Material =simulation{6}
Type_Material_Tube=simulation{7}
Sizing={A_cond, P_cond,[],[],[],Type_Material,Type_Material_Tube,Type_A,Type_B}

function log_mean = LMTD(T_CI,T_CO,T_HI,T_HO)
if (T_HI - T_CO) - (T_HO - T_CI) == 0;
    log_mean = abs(T_HI - T_CO);
else
    log_mean = ((T_HI - T_CO) - (T_HO - T_CI))/log((T_HI - T_CO)/(T_HO - T_CI));
end
end

end