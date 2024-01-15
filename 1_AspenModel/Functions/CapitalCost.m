function [CBM,CBM0]=CapitalCost(Sizing)

%%%%%%%%%%%%% EQUIPMENT COST CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   DATE: 09/JUL/2020                                                     %
%   AUTHOR: KELVIN PACHECO                                                %
%   EMAIL:  kelvinpac@gmail.com                                           %
%   VERSION: 1.0                                                          %
%   METHODOLOGY:                                                          %
%   The equipment module costing technique is a common technique to       %
%   estimate the cost of a new chemical plant (best for making preliminary%
%   cost estimates). It was introduced by Guthrie (1969) and forms the    %
%   basis of many of the equipment module techniques in use today.        %
%   The approach describes the costs of the purchase for the equipments on%
%   a base condition, a correction if different conditions (equipment type%
%   , pressure and construction materials) are used, a correction must be %
%   applied.                                                              %
%                                                                         %
%   INPUTS:                                                               %
%   Sizing = {Cap_Unit,   Actual value of the sizing [m2 or m3 or kW, etc]%
%             P,                Pressure [barg]                           %
%             D_Col,            Diameter of the column                    %
%             N_tray,           Number of trays                           %
%             N_Spares,         Number of spare equipments                %
%             Type_Material,    Material of construction                  %
%             Type_Material_Tube,   In the case of heat exchanger         %
%             Type_A,           Category of equipment (see 'equipments.mat%
%             Type_B}           Specific description of equipment         %
%                                                                         %
%   OUTPUT:                                                               %
%   [CBM,CBM0]=CapitalCost(Sizing)                                        %
%   where CBM is the Bare Module Cost (see Turton et al. 2018)            %
%         CBM0 is the bare module cost for ambient pressure               %
%   and carbon steel construction                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chemical Engineering Plant Cost Index (CEPCI) can be used to account for 
% changes that result from inflation.
% The overall CEPCI started 2019 at over 618 (618.7), and ended up the year
% under 600 (at 592.1), the yearly average for 2019 was 607.5.
% Source: www.chemengonline.com/2019-chemical-engineering-plant-cost-index-annual-average/

CEPCI=607.5;

% Retrieve the information from the cell called 'Sizing'

Cap_Unit            =double(Sizing{1});
P                   =Sizing{2};
D_Col               =Sizing{3};
N_tray              =Sizing{4};
N_Spares            =Sizing{5};
Type_Material       =Sizing{6};
Type_Material_Tube  =Sizing{7};
Type_A              =Sizing{8};
Type_B              =Sizing{9};

% Load the tables stored in Turton.mat (see Appendix A in Turton et al.
% 2018)

load('Turton.mat', 'TableA1', 'TableA2', 'TableA3_HE', 'TableA3_PV', 'TableA4', 'TableA6')

Unit_Row = (1:size(TableA1,1))*strcmp(TableA1(:,6),Type_B);
% *future implementation -> error when the capacity is exceeded.
% Cap_Min  = TableA1{Unit_Row,4};   
% Cap_Max  = TableA1{Unit_Row,5};
K1 = TableA1{Unit_Row,1};
K2 = TableA1{Unit_Row,2};
K3 = TableA1{Unit_Row,3};

% 1 . Purhased cost for base conditions
log_Cp0 = K1 + K2*log10(Cap_Unit) + K3*((log10(Cap_Unit))^2);
Cp0 = 10^log_Cp0;

% 2 . Update the cost to 2020
Cp0_currentyear = Cp0*(CEPCI/397);

% 3 . For exchangers, pumps, and vessels, find the pressure factor, 
%     Fp , Table A.2 and Equation (A.2) or (A.3)

if strcmp(Type_A,'Vessel')
    Fp = ( ((P+1)*D_Col/(2*(849.6-0.6*(P+1)))) + 0.00315 )/0.0063;
elseif strcmp(Type_A,'Heat Exchanger') || strcmp(Type_A,'Heater') || strcmp(Type_A,'Pump')
    P_min = TableA2{Unit_Row,4};
    if P > P_min
        C1 = TableA2{Unit_Row,1};
        C2 = TableA2{Unit_Row,2};
        C3 = TableA2{Unit_Row,3};
        log_Fp = C1 + C2*log10(P)+C3*log10(P)^2;
        Fp = 10^log_Fp;
    else
        Fp=1;
    end
    else
    Fp=1;
end

% 3a . For exchangers, pumps, and vessels, find the material of construction factor, 
% FM , Equation (A.4), Table A.3, and Figure A.18

if strcmp(Type_A,'Heat Exchanger')
    M_Row = (1:size(TableA3_HE,1))*strcmp(TableA3_HE(:,10),Type_B);
    M_Column= (1:size(TableA3_HE,2))*(strcmp(TableA3_HE(1,:),Type_Material) & strcmp(TableA3_HE(2,:),Type_Material_Tube))';
    FM = TableA3_HE{M_Row,M_Column};
    B1 = TableA4{Unit_Row,1};
    B2 = TableA4{Unit_Row,2};
    FBM = (B1+B2*FM*Fp);
    CBM = Cp0_currentyear*FBM;
    CBM0 = Cp0_currentyear*(B1+B2*1);
elseif strcmp(Type_A,'Pump') || strcmp(Type_A,'Vessel')
    M_Row    = (1:size(TableA3_PV,1))*strcmp(TableA3_PV(:,10),Type_B);
    M_Column = (1:size(TableA3_PV,2))*strcmp(TableA3_PV(1,:),Type_Material)';
    FM = TableA3_PV{M_Row,M_Column};
    B1 = TableA4{Unit_Row,1};
    B2 = TableA4{Unit_Row,2};
    FBM = (B1+B2*FM*Fp);
    if strcmp(Type_A,'Pump')
        CBM = Cp0_currentyear*(N_Spares+1)*FBM;
        CBM0 = Cp0_currentyear*(N_Spares+1)*(B1+B2*1);
    else
        CBM = Cp0_currentyear*FBM;
        CBM0 = Cp0_currentyear*(B1+B2*1);
        if (FM<1)||(Fp<1)
            CBM=CBM0;
        end
    end
end

% 3b . Find the correct relationship for bare module factor

if strcmp(Type_A,'Compressor') || strcmp(Type_A,'Turbine') || strcmp(Type_A,'Tray') 
    M_Row    = (1:size(TableA6,1))*strcmp(TableA6(:,4),Type_B);
    M_Column = (1:size(TableA6,2))*strcmp(TableA6(1,:),Type_Material)';
    FBM = TableA6{M_Row,M_Column};
    if strcmp(Type_A,'Tray')
        if (N_tray<20)
            log_Fq = 0.4771 + 0.08516*log10(N_tray)-0.3473*(log10(N_tray)^2);
            Fq=10^log_Fq;
        else
            Fq=1;
        end
        CBM = Cp0_currentyear*N_tray*FBM*Fq;
        CBM0=Cp0_currentyear*N_tray;
    else
        CBM = Cp0_currentyear*FBM;
        CBM0=Cp0_currentyear*FBM;
    end
elseif strcmp(Type_A,'Heater') || strcmp(Type_A,'Mixer') || strcmp(Type_A,'Reactor') || strcmp(Type_A,'Packing')
    FBM = TableA6{Unit_Row+1,1};
    CBM = Cp0_currentyear*N_tray*FBM;
    CBM0= Cp0_currentyear*3.24;
end

end