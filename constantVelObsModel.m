function [obsStateNext] = constantVelObsModel(obsStatePrev,dt)
    F_obs = [1, 0, dt, 0;
         0, 1, 0, dt;
         0, 0, 1, 0; 
         0, 0, 0, 1];
    obsStateNext = F_obs*obsStatePrev;
end