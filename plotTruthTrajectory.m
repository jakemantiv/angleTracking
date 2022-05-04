function [] = plotTruthTrajectory(truthStruct,ithMC)
if numel(truthStruct) == 1 % plot one trajectory at a time
xTgtTrue = truthStruct{1}.xTgtTrue;
xObsTrue = truthStruct{1}.xObsTrue;
figure(); hold on;
plot(xTgtTrue(1,:), xTgtTrue(2,:),'r');
plot(xObsTrue(1,:), xObsTrue(2,:),'b');
axis equal;
legend('Target Trajectory - Truth', 'Observer Trajectory - Truth', 'location', 'best');
title(['Engagement Ground Track - ', num2str(ithMC), ' MC']);
xlabel('East (m)'); ylabel('North (m)');

else % plot all on one plot
        figure(); hold on;
    for i = 1:numel(truthStruct)
    xTgtTrue = truthStruct{i}.xTgtTrue;
    xObsTrue = truthStruct{i}.xObsTrue;
    plot(xTgtTrue(1,:), xTgtTrue(2,:),'r');
    plot(xObsTrue(1,:), xObsTrue(2,:),'b');
    end
    axis equal;
    legend('Target Trajectory - Truth', 'Observer Trajectory - Truth', 'location', 'best');
    title(['Engagement Ground Track - ', num2str(numel(truthStruct)), ' MCs']);
    xlabel('East (m)'); ylabel('North (m)');
end 
end

