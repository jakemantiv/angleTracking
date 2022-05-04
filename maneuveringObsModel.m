function [obsStateNext, U] = maneuveringObsModel(obsStatePrev,dt, t, legTime, legHeading, legVel)

    R_turn = 150;
    
    notFoundLeg = true;
    legIdx = 0;
    for i = numel(legTime):-1:1
        if t > legTime(i) && notFoundLeg
            notFoundLeg = false;
            legIdx = i;
        end
    end
    currHeadingCmd = legHeading(legIdx);
    currVelCmd = legVel(:,legIdx);
    currLegTime = legTime(legIdx);
    timeSinceLastCmd = t - legTime(legIdx);
    M_time = R_turn/legVel(legIdx);
    turnRate = 360/M_time;
    prevHeading = atan2d(obsStatePrev(4), obsStatePrev(3));
    turnDone = abs(prevHeading - currHeadingCmd) < 1e-16;
    
    if ~turnDone
        headingErr = currHeadingCmd-prevHeading;
        if abs(headingErr) <= dt*turnRate
            newHeading = currHeadingCmd;
        else
            newHeading = prevHeading + dt*turnRate*sign(headingErr);
        end
        newVel = [cosd(newHeading)*currVelCmd; sind(newHeading)*currVelCmd];
        newPos = newVel*dt + obsStatePrev(1:2);
        obsStateNext = [newPos; newVel];
    else
         F_lin = [1, 0, dt, 0;
         0, 1, 0, dt;
         0, 0, 1, 0; 
         0, 0, 0, 1];
        obsStateNext = F_lin*obsStatePrev;
    end
    U = [obsStateNext(1) - obsStatePrev(1) - dt*obsStatePrev(3); 
         obsStateNext(2) - obsStatePrev(2) - dt*obsStatePrev(4);
         obsStateNext(3) - obsStatePrev(3);
         obsStateNext(4) - obsStatePrev(4)];
            
end