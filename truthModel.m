function [outStruct] = truthModel(truthModelInputStruct, target_mode)

Q = truthModelInputStruct.Q; 
R = truthModelInputStruct.R;
Gamma = truthModelInputStruct.Gamma;
xTrue0 = truthModelInputStruct.xTrue0; 
xTgt0 = truthModelInputStruct.xTgt0;
xObs0 = truthModelInputStruct.xObs0;
P0 = truthModelInputStruct.P0;
timeVec = truthModelInputStruct.timeVec;
legTime = truthModelInputStruct.obsLegTime;
legHeading = truthModelInputStruct.obsLegHeading;
legVel = truthModelInputStruct.obsLegVel;
dt = timeVec(2) - timeVec(1);
F = truthModelInputStruct.F;
xTrue = zeros(4,numel(timeVec));
xObsTrue = zeros(4,numel(timeVec));
xTgtTrue = zeros(4,numel(timeVec));
U = zeros(4,numel(timeVec));
z = zeros(numel(timeVec),1);
am = truthModelInputStruct.am;
xTrue(:,1) = xTrue0;
xObsTrue(:,1) = xObs0;
xTgtTrue(:,1) = xTgt0;
F_t = zeros(4,4,numel(timeVec));
Sw = chol(Q, 'lower'); 
Sv = chol(R,'lower');
Sp = chol(P0, 'lower');
for i = 1:numel(timeVec)
    if i == 1
        x0_draw = Sp*randn(4,1);
        prevXObsTrue = xObs0;
        xTruePrev = xTrue0 + x0_draw;
    else
        prevXObsTrue = xObsTrue(:,i-1);
        xTruePrev = xTrue(:,i-1);
    end
    
%     xObsTrue(:,i) = constantVelObsModel(prevXObsTrue,dt);
    [xObsTrue(:,i), U(:,i)] = maneuveringObsModel(prevXObsTrue,dt, timeVec(i),legTime,legHeading,legVel);
    
    x_draw = Sw*randn(2,1);
    
    if strcmpi(target_mode, 'straight')
        F_t(:,:,i) = F;
        xTrue(:,i) = F_t(:,:,i)*xTruePrev + Gamma*x_draw - U(:,i);
    elseif strcmpi(target_mode, 'clockwise')
        omega = -am/sqrt((xTruePrev(3) +prevXObsTrue(3))^2 + (xTruePrev(4) +prevXObsTrue(4))^2);
        F_t(:,:,i) = [1, 0, sin(omega*dt)/omega, -(1-cos(omega*dt))/omega;
               0, 1, (1-cos(omega*dt))/omega, sin(omega*dt)/omega;
               0, 0, cos(omega*dt), -sin(omega*dt);
               0, 0, sin(omega*dt), cos(omega*dt)];
        xTrue(:,i) = F_t(:,:,i)*(xTruePrev+prevXObsTrue) - xObsTrue(:,i) + Gamma*x_draw;
    elseif strcmpi(target_mode, 'counterclockwise')
        omega = am/sqrt((xTruePrev(3) +prevXObsTrue(3))^2 + (xTruePrev(4) +prevXObsTrue(4))^2);
        F_t(:,:,i) = [1, 0, sin(omega*dt)/omega, -(1-cos(omega*dt))/omega;
               0, 1, (1-cos(omega*dt))/omega, sin(omega*dt)/omega;
               0, 0, cos(omega*dt), -sin(omega*dt);
               0, 0, sin(omega*dt), cos(omega*dt)];
        xTrue(:,i) = F_t(:,:,i)*(xTruePrev+prevXObsTrue) - xObsTrue(:,i) + Gamma*x_draw;
    end
    
    xTgtTrue(:,i) = xTrue(:,i) + xObsTrue(:,i);
    
    z(i) = atan2(xTrue(2,i), xTrue(1,i)) + Sv*randn(1);
end
outStruct.z = z;
outStruct.xTrue = xTrue; 
outStruct.xTgtTrue = xTgtTrue;
outStruct.xObsTrue = xObsTrue;
outStruct.timeVec = timeVec;
outStruct.F_t = F_t;
outStruct.U = U;
end
