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

% Simulation Testing
Tsim = 120;
dt = 1;
Nsim = 100;

% favorable geometry
obsPos0 = [-300; 10]; 
obsHdg0 = 10;
obsVelMag0 = 15;
obsVel0 = [cosd(obsHdg0)*obsVelMag0; sind(obsHdg0)*obsVelMag0];

% unfavorable geometry
% obsPos0 = [-10; 10]; 
% obsVel0 = [1; 0];


tgtPos0 = [50; 90];
tgtVel0 = [1; 6];

xTgt0 = [tgtPos0; tgtVel0];
xObs0 = [obsPos0; obsVel0];

xTrue0 = xTgt0 - xObs0;

maxTime = 100;
obsLegTime = [0,25, 50, 75];
obsLegHeading = [obsHdg0,  165, 15, 165];
obsLegVel = [obsVelMag0, 15, 15, 15];

             

Gamma = [dt^2/2, 0;
         0,  dt^2/2;
         dt, 0
         0, dt];
F = [1, 0, dt, 0;
     0, 1, 0, dt;
     0, 0, 1, 0; 
     0, 0, 0, 1];
 
timeVec = dt:dt:maxTime;
N = numel(timeVec);

qW = .01;
Q = qW.*eye(2);
Q_PF = 0.015.*eye(2);
R = deg2rad(0.1);

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



constantVelPFinputStruct.R = R;
constantVelPFinputStruct.Q = Q_PF;
constantVelPFinputStruct.F = F;
constantVelPFinputStruct.Gamma = Gamma;
constantVelPFinputStruct.timeVec = timeVec;
constantVelPFinputStruct.Ns = 1000;
constantVelPFinputStruct.P0 = [100, 0, 0, 0; 
                              0, 100, 0, 0;
                              0, 0, 10, 0;
                              0, 0, 0, 10];
                          
constantVelPFinputStruct.xhat0 = xTrue0; % TODO: vary the initial state

for mc = 1:Nsim
    % Initial condition
%     delX0 = randn(n,1).*X_var;
%     X0 = delX0 + Xnom(0);
    
    % Truth model Simulation
    [truthStruct{mc}] = truthModel(truthModelInputStruct);
    
%     plotTruthTrajectory(truthStruct{mc}, mc);
    constantVelPFinputStruct.z = truthStruct{mc}.z;
    constantVelPFinputStruct.U = truthStruct{mc}.U;
    constantVelPFinputStruct.xObsTrue = truthStruct{mc}.xObsTrue;

    constantVelPFoutputStruct{mc} = constantVelPF(constantVelPFinputStruct);
    
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
      ex_PF(mc,:) = NEES(truthStruct{mc}.xTrue,constantVelPFoutputStruct{mc}.xhat_MMSE,constantVelPFoutputStruct{mc}.P);
%     ey_lkf(i,:) = NIS(Y,Yh_lkf,S_lkf);
%     ex_ekf(i,:) = NEES(X,Xh_ekf,P_ekf);
%     ey_ekf(i,:) = NIS(Y,Yh_ekf,S_ekf);
%     ex_ukf(i,:) = NEES(X,Xh_ukf,P_ukf);
%     ey_ukf(i,:) = NIS(Y,Yh_ukf,S_ukf);
    
    % Progress Notification
    fprintf('Done with %d\n',mc)
end
plotTruthTrajectory(truthStruct, mc)
analyzeOneRun(truthStruct{mc}, constantVelPFoutputStruct{mc});

alpha = 0.05;
exb_PF = mean(ex_PF,1);
n = 4;
rx = [chi2inv(alpha/2,Nsim*n), chi2inv(1 - alpha/2,Nsim*n)]'/Nsim;

% print some outputs of NEES/NIS tests to screen
% ratio of samples that were within ry bounds, should equal 1-alpha
nees_success_PF = sum(exb_PF<rx(2) & exb_PF>rx(1))./numel(exb_PF);


fprintf('PF NEES Test:\n');
fprintf('Observed ratio between bounds: %f\n', nees_success_PF);
fprintf('Expected ratio between bounds: %f\n\n', 1-alpha);

figure(); hold on;
plot(timeVec(2:end),exb_PF)
plot(timeVec, ones(size(timeVec)).*rx(1));
plot(timeVec, ones(size(timeVec)).*rx(2));




% NEES and NIS averaged over each time step
