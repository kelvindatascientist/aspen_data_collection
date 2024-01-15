function [N_Col, V_Col, N_Tray, A_Tray, D_Col]=Size_RadFrac_Vessel(simulation)

%%%%%%%%%%%%% EQUIPMENT SIZING: VAPOR-LIQUID CONTACT COLUMNS %%%%%%%%%%%%%%
%                                                                         %
%   OBJECTIVE:                                                            %
%                                                                         %
%   [Type, N_Col, V_Col, N_Tray, A_Tray, D_Col, L_Col, Errors]=           %
%   VLColumnSize(FV_V, MasD_V, NStage_Col) provides information           %
%   necessary to conduct cost calculations for vapor-liquid contact       %
%   columns.                                                              %
%                                                                         %
%   INPUT: The following values have to be provided:                      %
%                                                                         %
%       rho_L_top    : Top liquid stream mass density [kg/m^3]                    %
%       rho_V_top    : Top vapor stream mass density [kg/m^3]
%       V_g_top      : Top Vapor Rate [kg/s]
%       NStage_Col   : Number of equilibrium stages in the column [] 
%       rho_L_boilup : boilup liquid stream mass density [kg/m^3]                    %
%       rho_V_boilup : boilup vapor stream mass density [kg/m^3]
%       V_g_boilup   : boilup Vapor Rate [kg/s]       
%                                                                         %
%   OUTPUT: The function generates the following values:                  %
%                                                                         %
%       Type    : Equipment type:                                         %
%                 'V1' = Horizontal vessel                                %
%                 'V2' = Vertical vessel                                  %
%       N_Col   : Totla number of columns (series and parallel) used to   %
%                 satisfy the separation requirements.                    %
%       V_Col   : Column volume [m^3]                                     %
%       N_Tray  : Number of column trays                                  %
%       D_Col   : Column diameter [m]                                     %
%       L_Col   : Column length [m]                                       %
%       Errors  : Contains list of error identified during execution      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho_L_top    = simulation{1};
rho_V_top    = simulation{2};
V_g_top      = simulation{3}/3600;
rho_L_boilup = simulation{4};
rho_V_boilup = simulation{5};
V_g_boilup   = simulation{6}/3600;
N_Stages     = simulation{7};
%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT ERROR CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%
    Errors ={};
    if (V_g_boilup<0); Errors = [Errors; 'Vapor volume flow must be >=0']; end
    if (rho_V_boilup<=0); Errors =[Errors; 'Vapor mass density must be >0']; end
    if (N_Stages<2 || mod(N_Stages,1)~=0);Errors =...
    [Errors; 'Number of stages should be integer >=2']; end
    if size(Errors,1)>0;
        N_Col=0; V_Col=0; N_Tray=0;
        A_Tray=0; D_Col=0; L_Col=0; return;
    end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S = 0.6096; % or 24''  plate spacing between the trays  [m]  \cite{Biegler1999} pg.118                   %
    eff = 0.8;  %  plate efficiency [] \cite{Biegler1999} pg 121
    % The maximum column diameter should be 20 ft or 6 m \cite{Biegler1999} pg
    % 121
    % for calculations of cost using Turton required the maximum area of
    % 12.3m²
    D_max = 5;
    L_Max = 60;                 % maximum column length [m] \cite{Barker2018}
    LDr_Max = 30;               % Maximum L/D ratio	\cite{Barker2018}
    V_Min = 0.3;                % Minimum column volume \cite{Turton2012} annex A
    V_Max  = 520;               % Maximum column volume \cite{Turton2012} annex A
%%%%%%%%%%%%%%%%%%%%%%%% CALCULATIONS AND RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%

% Calculations for the tower
    u_w_top =    (-0.171*(S^2)+0.27*S-0.047)*(((rho_L_top-rho_V_top)/rho_V_top)^(1/2));
    D_top =    ((4*V_g_top)/(pi*rho_V_top*u_w_top))^(1/2);
    u_w_boilup = (-0.171*(S^2)+0.27*S-0.047)*(((rho_L_boilup-rho_V_boilup)/rho_V_boilup)^(1/2));
    D_boilup = ((4*V_g_boilup)/(pi*rho_V_boilup*u_w_boilup))^(1/2);
    
    D_Col = max(D_top,D_boilup);
%Calculation for trays
    A_Req = pi*(D_Col^2)/4;
    N_Tray = (N_Stages) ;      % Number of column trays without condenser and reboiler
    L_Req = (N_Stages/eff) * S;
    
    NSeries_Col = ceil(L_Req/L_Max) ;   % Number of columns in series
    
    while true
        L_Col = L_Req/NSeries_Col   ;   % Column length
        A_Max = V_Max/L_Col;           % Max. column cross area
        A_Min = V_Min/L_Col;            % Min. column cross area
        NParallel_Col = ceil(A_Req/A_Max) ; % Number of columns in parallel
        %NParallel_Col = ceil(D_Col/sqrt(4*12.3/pi))
        A_Col = max(A_Req/NParallel_Col,A_Min); % Column area
        D_Col = sqrt(4*A_Col/pi)   ;      % Column diameter
        if (L_Col/D_Col)>LDr_Max
            NSeries_Col = NSeries_Col+1;
        else
            break;
        end
    end
N_Col = NSeries_Col*NParallel_Col;  % Total number of columns
V_Col = A_Col*L_Col;                % Column volume
A_Tray = double(A_Col);

%second implementation where a maximum area is used only for trays
Total_Area = A_Tray*N_Tray;
if (Total_Area) > (12.3*N_Tray)
    A_Tray = 12.3;
    N_Tray = Total_Area/12.3;
end
    
        
end