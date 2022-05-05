function [F_t] = calcF(inputStruct, prevXObsTrue,xPrev ,target_mode)
    F = inputStruct.F;
    am = inputStruct.am;
    dt = inputStruct.timeVec(2) - inputStruct.timeVec(1);
    if strcmpi(target_mode, 'straight')
        F_t = F;
    elseif strcmpi(target_mode, 'clockwise')
        omega = -am/sqrt((xPrev(3) +prevXObsTrue(3))^2 + (xPrev(4) +prevXObsTrue(4))^2);
        F_t = [1, 0, sin(omega*dt)/omega, -(1-cos(omega*dt))/omega;
               0, 1, (1-cos(omega*dt))/omega, sin(omega*dt)/omega;
               0, 0, cos(omega*dt), -sin(omega*dt);
               0, 0, sin(omega*dt), cos(omega*dt)];
    elseif strcmpi(target_mode, 'counterclockwise')
        omega = am/sqrt((xPrev(3) +prevXObsTrue(3))^2 + (xPrev(4) +prevXObsTrue(4))^2);
        F_t = [1, 0, sin(omega*dt)/omega, -(1-cos(omega*dt))/omega;
               0, 1, (1-cos(omega*dt))/omega, sin(omega*dt)/omega;
               0, 0, cos(omega*dt), -sin(omega*dt);
               0, 0, sin(omega*dt), cos(omega*dt)];
    end
end

