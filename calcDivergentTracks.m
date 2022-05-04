function [divergent_idx] = calcDivergentTracks(truthStruct, constantVelPFoutputStruct, posErrThreshold, timeDivergedThreshold)

divergent_idx = [];
for mc = 1:numel(truthStruct)
    startTimeDiverged = -1;
    timeDiverged = 0;
    divergenceConfirmed = false;
    for i = 1:numel(truthStruct{mc}.timeVec)
        estErr = truthStruct{mc}.xTrue(:,i) - constantVelPFoutputStruct{mc}.xhat_MMSE(:,i);
        posErr = sqrt(estErr(1)^2 + estErr(2)^2);
        if posErr >= posErrThreshold
            if startTimeDiverged == -1
                timeDiverged = 0;
                startTimeDiverged = truthStruct{mc}.timeVec(i);
            else
                timeDiverged = truthStruct{mc}.timeVec(i) - startTimeDiverged;
            end
        else
            startTimeDiverged = -1;
        end
        if timeDiverged >= timeDivergedThreshold && ~divergenceConfirmed
           divergent_idx = [divergent_idx; mc]; 
           divergenceConfirmed = true;
        end
    end
    
end
end