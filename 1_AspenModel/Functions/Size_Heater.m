function [Sizing,N_HeatEx,A_req]=Size_Heater(simulation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_in   = simulation{2};
T_out  = simulation{3};
Q_duty = simulation{4}; % W
Press  = simulation{5}; % bar
if strcmp(simulation{1},'Cooler')
    Type_A = 'Heat Exchanger';
    Type_B = 'Heat Exchanger U tube';
    % o condensador é um trocador de casco e tubo feito de aço carbono
    % refrigeração água a 4bar e 30ºC, saindo do trocador a 40ºC
    U = 765;  %W/m2.K, tabela do Campbell (2011) Figure 13.2
    if T_out > 38
        % usar Cooling water
        T_supply = 25; %oC, considerando 4 bar
        T_return = 30; %oC , pode ser 50oC precisa de ref.
    elseif (T_out <= 38) && (T_out > -20)
        T_supply = -25; %oC
        T_return = -25; %oC
    else
        % usar Freon 11
        T_supply = -40; %oC
        T_return = -40; %oC 
    end
    A_req = -Q_duty/(U*LMTD(T_supply,T_return,T_in,T_out));
else
    % starts if is a heater
    Type_A = 'Heat Exchanger';
    Type_B = 'Heat Exchanger U tube';
    % refervedor é um kettle reboiler de aço carbono
    % aquecimento no refeverdor, utiliza-se, quando possível, vapor saturado de
    % baixa pressão - 6 bar (160ºC), ou de média pressão - 10 bar (184ºC)
    U = 910;  %W/m2.K, tabela do Campbell (2011) Figure 13.2
    T_supply = 160; %oC, considerando vapor fraction 1
    T_return = 160; %oC , considerando vapor fraction 0
    LMTD(T_in,T_out,T_supply,T_return);
    A_req = Q_duty/(U*LMTD(T_in,T_out,T_supply,T_return));
end

A_min  = 10;
A_max  = 1000;

%%% Calculating the minimum number of units
N_HeatEx = ceil(A_req/A_max);

%%% Sizing
A_Prelim = A_req/N_HeatEx;
A = max(A_Prelim,A_min);
% Sizing = {Cap_Unit, P, D_Col,N_tray, N_Spares,Type_Material,Type_Material_Tube,Type_A,Type_B}
Type_Material =simulation{6};
Type_Material_Tube=simulation{7};
Sizing={A, Press,[],[],[],Type_Material,Type_Material_Tube,Type_A,Type_B};
N_HeatEx = ceil(A_req/A_max);
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