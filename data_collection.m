%% Example of code for data collection
%
% Kelvin Pacheco
%
% Phd Student - University of Sao Paulo - Department of Chemical
% Engineering
%
% Date: 11/05/20
%% DATA RETRIEVE
% The files that must be in the folder:
% Simulation of Aspen (Simulation_Name)
clc
clear all

% The first step is initialize the application, controlling Aspen Plus using ActiveX in MATLAB
% Link original do https://www.mathworks.com/matlabcentral/fileexchange/69464-aspen-plus-matlab-link?focused=1eefd5fe-a0ce-4f86-8f7e-debdcabd5760&tab=function

Aspen = actxserver('Apwn.Document.38.0');     % Initialize suit (34.0 -> V8.8; 35.0 -> V9.0; 36.0 -> V10.0)
[stat,mess]=fileattrib;                       % get attributes of folder (Necessary to establish the location of the simulation)
Simulation_Name = 'CabaleroV12';              % Aspeen Plus Simulation Name
Aspen.invoke('InitFromArchive2',[mess.Name '\1_AspenModel\Simulation\' Simulation_Name '.bkp']);
Aspen.Visible = 1;                            % 1 -> Aspen is Visible; 0 -> Aspen is open but not visible
Aspen.SuppressDialogs = 1;                    % Suppress windows dialogs.

% Run the Aspen Plus Simulation

Aspen.Engine.Run2(1);
while Aspen.Engine.IsRunning == 1             % 1 -> If Aspen is running; 0 -> If Aspen stop.
    pause(0.5);
end

%% Application in real problems

% Define the upper and lower bound of the sampling
lb  =  [1.5; 0.15];
ub  =  [2.5; 0.25];

% Use the latin hypercube to sampling
design = make_design(lb,ub,2,'lhsd')
%%
% Data collection for a RadFrac
% Insert the name of a related block in Aspen Plus, feed stream and output
% streams
Block_Name='B1';
Feed_Stream='101';
Distillate_Stream='S4';
Bottom_Stream='101';

% Run a 'for' loop taking into account the DoE previously performed. In the
% case of a variable is needed this site
% (https://chejunkie.com/knowledge-base/navigating-variable-explorer-aspen-plus/)
% provides the guidelines

for i = 1:size(design,2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assign design of experiment values to specific variables in Aspen   %   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Input\BASIS_RR']).Value        = design(1,i);   % Column Reflux
    Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Input\D:F']).Value             = design(2,i);            % Distillate to feed ratio
    %Aspen.Tree.FindNode("\Data\Blocks\B1\Input\NSTAGE").Value          = fix(design(3,i));       % Number of stages
    %Aspen.Tree.FindNode("\Data\Blocks\B1\Input\FEED_STAGE\S1").Value = fix(design(3,i)*design(4,i));     % Stage feed
    %Aspen.Tree.FindNode("\Data\Blocks\CL01\Input\PRES1").Value           = Column_pressure(i);% Column Pressure
    %Aspen.Tree.FindNode("\Data\Streams\FEED\Input\TEMP").Value           = Temp_feed(i);  
    %Aspen.Tree.FindNode("\Data\Streams\FEED\Input\FLOW\MIXED\H2").Value      = Z_a;
    %Aspen.Tree.FindNode("\Data\Streams\FEED\Input\FLOW\MIXED\METHANOL").Value= Z_b;
    %Aspen.Tree.FindNode("\Data\Streams\FEED\Input\FLOW\MIXED\ACETIC").Value  = Z_c;
    %Aspen.Tree.FindNode("\Data\Streams\FEED\Input\FLOW\MIXED\CO2").Value     = Z_d;
    %Aspen.Tree.FindNode("\Data\Streams\FEED\Input\FLOW\MIXED\WATER").Value   = Z_e;
    %Aspen.Tree.FindNode("\Data\Streams\FEED\Input\FLOW\MIXED\CH4").Value     = Z_f;
    %Aspen.Tree.FindNode("\Data\Streams\FEED\Input\FLOW\MIXED\METHY-01").Value= Z_g;
    %Aspen.Tree.FindNode("\Data\Streams\FEED\Input\FLOW\MIXED\ETHYL-01").Value= Z_h;
    %Aspen.Tree.FindNode("\Data\Streams\FEED\Input\FLOW\MIXED\ETHAN-01").Value= Z_i;
    %Aspen.Tree.FindNode("\Data\Streams\FEED\Input\FLOW\MIXED\SOLVENT").Value = Z_j;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reinitialize the simulation                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Aspen.Reinit;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run the simulation                                                  %
    % (1) ---> Matlab isn't busy                                          %
    % (0) ---> Matlab is busy                                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Aspen.Engine.Run2(1);
    time = 1;
    while Aspen.Engine.IsRunning == 1 
        % 1-> If Aspen is running; 0-> If Aspen stop.
        pause(0.5);
        time = time+1;
        if time==350 
            %Control of simulation time.
            Aspen.Engine.Stop;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % End of the simulation, check the convergence                        %
    % (1) ---> Aspen Plus didn't converge                                 %
    % (0) ---> Aspen Plus converged                                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Simulation_Convergency = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Output\BLKSTAT']).Value    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save the results                                                    %
    % Duty_reb : Duty heat in reboiler [W]
    % simulation{1}=rho_L_top : Top liquid stream mass density [kg/m^3]
    % simulation{2}=rho_V_top : Top vapor stream mass density [kg/m^3]
    % simulation{3}=V_g_top   : Top Vapor Rate [kg/s]
    % simulation{4}=NStage_Col   : Number of equilibrium stages in the column [] 
    % simulation{5}=rho_L_boilup : boilup liquid stream mass density [kg/m^3]    
    % simulation{6}=rho_V_boilup : boilup vapor stream mass density [kg/m^3]
    % simulation{7}=V_g_boilup   : boilup Vapor Rate [kg/s]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    if Simulation_Convergency == 0 && time < 1000
        Duty_reb(i)   = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Output\REB_DUTY']).Value; 
        Duty_cond(i)  = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Output\COND_DUTY']).Value;
        Pent(i)       = Aspen.Tree.FindNode(['\Data\Streams\',Distillate_Stream,'\Output\MOLEFRAC\MIXED\PENTANE']).Value
        MoleFlow(i)   = Aspen.Tree.FindNode(['\Data\Streams\',Distillate_Stream,'\Output\MOLEFLOW\MIXED\PENTANE']).Value
        N_Stage(i)    = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Input\NSTAGE']).Value
        simulation{1} = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Output\HYD_RHOL\1']).Value;
        simulation{2} = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Output\HYD_RHOV\1']).Value;
        simulation{3} = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Output\HYD_VMF\1']).Value;
        simulation{4} = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Output\HYD_RHOL\',num2str(N_Stage(i)-1)]).Value;
        simulation{5} = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Output\HYD_RHOV\',num2str(N_Stage(i)-1)]).Value;
        simulation{6} = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Output\HYD_VMF\',num2str(N_Stage(i)-1)]).Value;
        simulation{7} = double(N_Stage(i));
        % Sizing the vessel
        cd ([mess.Name '\1_AspenModel\Functions\'])
        [N_Col, V_Col, N_Tray, A_Tray, D_Col]=Size_RadFrac_Vessel(simulation)
        P_tower= Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Input\PRES1']).Value
        
        % Capital Cost for vessel
        load('equipments.mat')
        Number_corresponding=35;
        Sizing = {V_Col, P_tower, D_Col,N_Tray,[],'CS','CS',equipments{Number_corresponding,1,1},equipments{Number_corresponding,2,1}}
        [A,A0]=CapitalCost(Sizing);
        CBM_Vess(i)=A;
        CBM0_Vess(i)=A0;
        
        % Capital Cost for trays
        Number_corresponding=37;
        Sizing = {A_Tray, P_tower, D_Col,N_Tray,[],'CS','CS',equipments{Number_corresponding,1,1},equipments{Number_corresponding,2,1}}
        [B,B0]=CapitalCost(Sizing);
        CBM_Tray(i)=B;
        CBM0_Tray(i)=B0;
        cd ([mess.Name])
        
        % Collect data to size and costing the condenser and reboiler
        simulationB{1} = Duty_cond(i)/N_Col;
        simulationB{2} = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Output\B_TEMP\1']).Value;%T_top
        simulationB{3} = P_tower; %P_cond
        simulationB{4} = Duty_reb(i)/N_Col;
        simulationB{5} = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Output\HYD_TL\',num2str(N_Stage(i))]).Value; %T_bot
        simulationB{6} = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Output\HYD_TL\',num2str(N_Stage(i)-1)]).Value;
        simulationB{7} = Aspen.Tree.FindNode(['\Data\Blocks\',Block_Name,'\Output\B_PRES\',num2str(N_Stage(i))]).Value;%P_bot
        cd ([mess.Name '\1_AspenModel\Functions\'])
        [C,C0,D,D0]=Size_RadFRac_condreb(simulationB);
        cd ([mess.Name])
        CBM_cond(i) = C;
        CBM0_cond(i)= C0;
        CBM_reb(i)  = D;
        CBM0_reb(i) = D0;
        
        % Calculate the total cost of the tower (Vessel+Tray+Cond+Reb)
        CBM_Tower(i) = CBM_Vess(i)*N_Col + CBM_Tray(i)*N_Col + CBM_cond(i)*N_Col+CBM_reb(i)*N_Col
        CBM0_Tower(i)= CBM0_Vess(i)*N_Col+CBM_Tray(i)*N_Col+CBM0_cond(i)*N_Col+CBM0_reb(i)*N_Col
        
        % Retrieve information about the utilities
        % Water usage in condenser in [kg/h]
        % LP steam usage in reboiler in [kg/h]
        water_usage(i)=Aspen.Tree.FindNode(['\Data\Utilities\WATER\Output\UTL_USAGE\',Block_Name]).Value
        lp_usage(i)   =Aspen.Tree.FindNode(['\Data\Utilities\LP\Output\UTL_USAGE\',Block_Name]).Value
    
        % Clear some variables
        clear simulation
        clear simulationB
    else
        % This part represents the non convergence of the simulation,
        % therefore, all variables are assigned as zero
        Duty_reb(i)  = 0;
        Duty_cond(i) = 0;
        Pent(i)      = 0;
        MoleFlow(i)  = 0;
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close the Aspen Plus application                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Aspen.Close;
Aspen.Quit;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save the data collected in a dataset                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataset = design(:,1:size(Pent,2));
dataset(3,:) = Pent;
dataset(4,:) = Duty_reb;
dataset(5,:) = Duty_cond;
dataset(6,:) = MoleFlow;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean the dataset                                                       %
% During the retrieve of simulation, when there was a non convergence     %
% problem it was flaged as zero for all the results. In this case, these  %
% flaged rows will be removed.                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataset(:,(dataset(4,:)==0)) = [];

% For the neural networks each column represents a feature and each row a
% sample, therefore, the variable is transposed
dataset=dataset'
%%
% Save the files related into the Folder
NewFolder=[mess.Name '\2_Datasets\' Block_Name ];
if ~exist(NewFolder, 'dir')
    mkdir(NewFolder)
end

if numel(dir(NewFolder)) <=2
    fprintf('Empty Folder. \n')
    counter=1
else
    fprintf('NOT Empty Folder. \n')
end

if isfile([NewFolder,'\','dataset',num2str(counter),'_',Block_Name,'.mat'])
     fprintf('File already exists. \n')% File exists.
     counter=counter+1;
     save([NewFolder,'\','dataset',num2str(counter),'_',Block_Name,'.mat'],'dataset');
     csvwrite([NewFolder,'\','dataset',num2str(counter),'_',Block_Name,'.csv'],dataset);
else
     fprintf('File does not exist yet. \n')
     save([NewFolder,'\','dataset',num2str(counter),'_',Block_Name,'.mat'],'dataset');
     csvwrite([NewFolder,'\','dataset',num2str(counter),'_',Block_Name,'.csv'],dataset);
end

%%
input_test=dataset(1:2,:);
output_test=dataset(3:6,:);
csvwrite('input_test.csv', input_test);
csvwrite('output_test.csv', output_test);
%% merge datasets

merged_dataset=[dataset_AA;dataset_AB]
csvwrite('FileName.csv', merged_dataset);
%% Plotting
figure()

figure(1)
plot(design(1,:), Duty_reb, 'ok');
%plot(design(1,1:size(Pent,2)), Pent, 'ok');
xlabel('Reflux Ratio')
ylabel('Duty')

%%
%%%%  Auxiliary functions ###################################################
function [design] = make_design(lb,ub,nex,tipo)
% This function obtain the sample for input of the surrogate modelling
% inputs: lb, ub, nex and tipo

switch tipo
    case 'rand'
        norm = rand(nex, length(lb));
        design = norm'.*(ub-lb)+lb;
    case 'lhsd'
        norm = lhsdesign(nex, length(lb),'iteration',50); %https://www.mathworks.com/help/stats/lhsdesign.html
        design = norm'.*(ub-lb)+lb;
end
end