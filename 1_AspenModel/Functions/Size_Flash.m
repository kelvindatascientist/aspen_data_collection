function [V_Col,H_vess,D_vess,N_vess,Type_A,Type_B]=Size_Flash(simulation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% retrieve data from cell called 'simulation'
%simulation = {rho_V,rho_L,Q_v,Pressure,Number_corresponding}
load('equipments.mat')
rho_V                = simulation{1}; %kg/m3
rho_L                = simulation{2}; %kg/m3
Q_v                  = simulation{3}; %m3/s
Pressure             = simulation{4}; %bar
Number_corresponding = simulation{5}; %[]

% retrieve the name of the equipment
% Type_A = 'Vessel'
Type_A=equipments{Number_corresponding,1,1};
% Type_B = 'Vessel Vertical'
Type_B=equipments{Number_corresponding,2,1};

H_D = 3;

% retrieve the min and max capacity
load('Turton.mat')
Unit_Row = (1:size(TableA1,1))*strcmp(TableA1(:,6),Type_B);
Cap_Min  = TableA1{Unit_Row,4};
Cap_Max  = TableA1{Unit_Row,5};

K_s = 0.107;

% Eq. Souders Brown
vv = K_s*sqrt((rho_L - rho_V)/rho_V);

% eq. campbell (based on vapor velocity)
D_vess = sqrt(4*Q_v/(pi*rho_V*vv));
H_vess = D_vess*H_D;

A_Req = pi*(D_vess^2)/4;
V_Req = A_Req*H_vess;
N_vess = ceil(V_Req/Cap_Max);

A_vess = A_Req/N_vess;
D_vess = sqrt(4*A_vess/pi);
V_Col = V_Req/N_vess;
H_vess = V_Col/A_vess;
end