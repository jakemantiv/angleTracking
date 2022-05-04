clear; close all;
% clc;

% save figs? 
saveFigs = false;

% path to saved figures
figPath = ['..' filesep, 'figs', filesep];

% fix random seed
rng(12)

% where to save the data? 
dataFileName = 'mcTestData';

posErrThreshold = 150;
timeDivergedThreshold = 15;

% Simulation Testing
Tsim = 120;
dt = 0.5;
Nsim = 100;

% favorable geometry
obsPos0 = [-200; 0]; 
obsHdg0 = 20;
obsVelMag0 = 20; % favorable geometry

% obsVelMag0 = 0; % unfavorable geometry
 
obsVel0 = [cosd(obsHdg0)*obsVelMag0; sind(obsHdg0)*obsVelMag0];



tgtPos0 = [0; 100];
tgtVel0 = [0; 15];

xTgt0 = [tgtPos0; tgtVel0];
xObs0 = [obsPos0; obsVel0];

xTrue0 = xTgt0 - xObs0;

maxTime = 100;
obsLegTime = [0,25, 50, 75];
obsLegHeading = [obsHdg0,  160, 20, 160];
obsLegVel = [obsVelMag0, obsVelMag0, obsVelMag0, obsVelMag0];

             

Gamma = [dt^2/2, 0;
         0,  dt^2/2;
         dt, 0
         0, dt];
F = [1, 0, dt, 0;
     0, 1, 0, dt;
     0, 0, 1, 0; 
     0, 0, 0, 1];
P0 = [50, 0, 0, 0; 
      0, 50, 0, 0;
      0, 0, 2.5, 0;
      0, 0, 0, 2.5];
 
timeVec = dt:dt:maxTime;
N = numel(timeVec);

qW = .01;
Q = qW.*eye(2);
Q_PF = 0.010.*eye(2);
R = deg2rad(1);

xTrue = zeros(4,numel(timeVec));
xObsTrue = zeros(4,numel(timeVec));
xTgtTrue = zeros(4,numel(timeVec));
z = zeros(numel(timeVec),1);

% MonteCarlo Simulation Runs
ex_lkf = zeros(Nsim,N-1);
ey_lkf = zeros(Nsim,N-1);
ex_ekf = zeros(Nsim,N-1);
ey_ekf = zeros(Nsim,N-1);
ex_ukf = zeros(Nsim,N-1);
ey_ukf = zeros(Nsim,N-1);

truthModelInputStruct = struct();
truthModelInputStruct.Q = Q; 
truthModelInputStruct.R = R;
truthModelInputStruct.F = F;



truthModelInputStruct.Gamma = Gamma;
truthModelInputStruct.xTrue0 = xTrue0; 
truthModelInputStruct.xTgt0 = xTgt0;
truthModelInputStruct.xObs0 = xObs0;
truthModelInputStruct.timeVec = timeVec;
truthModelInputStruct.obsLegTime = obsLegTime;
truthModelInputStruct.obsLegHeading = obsLegHeading;
truthModelInputStruct.obsLegVel = obsLegVel;
truthModelInputStruct.am = 0.2; % lateral accel capability of target;
truthModelInputStruct.P0 = P0;

constantVelPFinputStruct.R = R;
constantVelPFinputStruct.Q = Q_PF;
constantVelPFinputStruct.F = F;
constantVelPFinputStruct.xObs0 = xObs0;
constantVelPFinputStruct.Gamma = Gamma;
constantVelPFinputStruct.timeVec = timeVec;
constantVelPFinputStruct.Ns = 1000;

constantVelPFinputStruct.P0 = 2.*P0;
                          
constantVelPFinputStruct.xhat0 = xTrue0; % TODO: vary the initial state
target_mode = 'clockwise'; % 'straight', 'clockwise', 'counterclockwise';
for mc = 1:Nsim
    % Initial condition
%     delX0 = randn(n,1).*X_var;
%     X0 = delX0 + Xnom(0);
    
    % Truth model Simulation
    [truthStruct{mc}] = truthModel(truthModelInputStruct,target_mode);
    
%     plotTruthTrajectory(truthStruct{mc}, mc);
    constantVelPFinputStruct.z = truthStruct{mc}.z;
    constantVelPFinputStruct.U = truthStruct{mc}.U;
    constantVelPFinputStruct.xObsTrue = truthStruct{mc}.xObsTrue;
    constantVelPFinputStruct.F_t = truthStruct{mc}.F_t;

    constantVelPFoutputStruct{mc} = singleTargetPF(constantVelPFinputStruct,target_mode);
    
%     constantVelPFoutputStruct{mc} = singleTargetPF(constantVelPFinputStruct,'clockwise');

    % Linearized Kalman Filter Estimate
%     [Xh_lkf,Yh_lkf,P_lkf,S_lkf,Sx_lkf] = KF(time,Y,U,X0,P0,Xnom,Unom,Ynom,F,G,H,M,Om,Q,R);
%     
%     % Extended Kalman Filter
%     [Xh_ekf,Yh_ekf,P_ekf,S_ekf,Sx_ekf] = EKF(time,Y,U,X0,P0,Xnom,Unom,Fnl,F,A,Hnl,H,Om,Q,R);
%     
%     % Unscented Kalman Filter
%     [Xh_ukf,Yh_ukf,P_ukf,S_ukf,Sx_ukf] = UKF(time,Y,U,X0,0.1.*P0,Fnl,Hnl,H,Om,Q,R);
% 
%     % Calculate NEES and NIS
%     ey_lkf(i,:) = NIS(Y,Yh_lkf,S_lkf);
%     ex_ekf(i,:) = NEES(X,Xh_ekf,P_ekf);
%     ey_ekf(i,:) = NIS(Y,Yh_ekf,S_ekf);
%     ex_ukf(i,:) = NEES(X,Xh_ukf,P_ukf);
%     ey_ukf(i,:) = NIS(Y,Yh_ukf,S_ukf);
    
    % Progress Notification
    fprintf('Done with %d\n',mc)
end
[divergent_idx] = calcDivergentTracks(truthStruct, constantVelPFoutputStruct, posErrThreshold, timeDivergedThreshold)
all_idx = (1:Nsim);
non_divergent_idx = all_idx(~ismember(all_idx,divergent_idx));
disp(['Num divergent tracks:', num2str(numel(divergent_idx))]); 
    
RMS = calcRMS(non_divergent_idx, truthStruct, constantVelPFoutputStruct);
% calcCRLBRMS
RTAMS = calcRTAMS(non_divergent_idx, truthStruct, constantVelPFoutputStruct, obsLegTime)

plotTruthTrajectory(truthStruct, mc)
analyzeOneRun(truthStruct{mc}, constantVelPFoutputStruct{mc});

figure(); 
plot(timeVec, RMS);




% NEES and NIS averaged over each time step
