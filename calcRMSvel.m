function [RMS] = calcRMSvel(non_divergent_idx, truthStructs, outputStructs)

M = numel(non_divergent_idx);
RMS = zeros(size(truthStructs{1}.xTrue,2),1);
for i = 1:numel(RMS)
    mySum = 0;
    for mc_idx = 1:numel(non_divergent_idx)
        mc = non_divergent_idx(mc_idx);
        estErr = truthStructs{mc}.xTrue(:,i) - outputStructs{mc}.xhat_MMSE(:,i);
        mySum = mySum + estErr(3)^2 + estErr(4)^2;
    end
    RMS(i) = sqrt((1/M)*mySum);
end

end