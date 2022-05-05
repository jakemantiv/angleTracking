function [RTAMS] = calcRTAMSvel(non_divergent_idx, truthStructs, outputStructs, obsLegTime)
M = numel(non_divergent_idx);
RTAMS= zeros(size(truthStructs{1}.xTrue,2),1);
tmax_idx = size(truthStructs{1}.xTrue,2);
tmax = truthStructs{1}.timeVec(tmax_idx);
l = obsLegTime(2);
l_idx = find(truthStructs{1}.timeVec >= l, 1,'first');
for i = (l_idx+1):tmax_idx
    mySum = 0;
    for mc_idx = 1:M
        mc = non_divergent_idx(mc_idx);
        estErr = truthStructs{mc}.xTrue(:,i) - outputStructs{mc}.xhat_MMSE(:,i);
        mySum = mySum + estErr(3)^2 + estErr(4)^2;
    end
    RTAMS(i) = sqrt((1/M)*mySum);
end

RTAMS = sqrt(1/(M*(tmax-l))*sum(RTAMS));
end