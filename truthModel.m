function [outStruct] = truthModel(truthModelInputStruct)

Q = truthModelInputStruct.Q; 
R = truthModelInputStruct.R;
F = truthModelInputStruct.F;
Gamma = truthModelInputStruct.Gamma;
xTrue0 = truthModelInputStruct.xTrue0; 
xTgt0 = truthModelInputStruct.xTgt0;
xObs0 = truthModelInputStruct.xObs0;
timeVec = truthModelInputStruct.timeVec;

dt = timeVec(2) - timeVec(1);

xTrue = zeros(4,numel(timeVec));
xObsTrue = zeros(4,numel(timeVec));
xTgtTrue = zeros(4,numel(timeVec));
z = zeros(numel(timeVec),1);

xTrue(:,1) = xTrue0;
xObsTrue(:,1) = xObs0;
xTgtTrue(:,1) = xTgt0;

Sw = chol(Q, 'lower'); 
Sv = chol(R,'lower');

for i = 1:numel(timeVec)
    if i == 1
        prevXObsTrue = xObs0;
        xTruePrev = xTrue0;
    else
        prevXObsTrue = xObsTrue(:,i-1);
        xTruePrev = xTrue(:,i-1);
    end
    
    xObsTrue(:,i) = constantVelObsModel(prevXObsTrue,dt);
    
    x_draw = Sw*randn(2,1);
    xTrue(:,i) = F*xTruePrev + Gamma*x_draw;
    
    xTgtTrue(:,i) = xTrue(:,i) + xObsTrue(:,i);
    
    z(i) = atan2(xTrue(2,i), xTrue(1,i)) + Sv*randn(1);
end
outStruct.z = z;
outStruct.xTrue = xTrue; 
outStruct.xTgtTrue = xTgtTrue;
outStruct.xObsTrue = xObsTrue;
outStruct.timeVec = timeVec;
end
