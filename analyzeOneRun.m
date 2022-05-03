function [] = analyzeOneRun(truthStruct, constantVelPFoutputStruct)

timeVec = truthStruct.timeVec;
z = truthStruct.z;
xTrue = truthStruct.xTrue;
xObsTrue = truthStruct.xObsTrue;
xTgtTrue = truthStruct.xTgtTrue;
xhat_tgt_MMSE =constantVelPFoutputStruct.xhat_tgt_MMSE;
w = constantVelPFoutputStruct.w;
xsamps_post = constantVelPFoutputStruct.xsamps_post;
xhat_MMSE = constantVelPFoutputStruct.xhat_MMSE;
P = constantVelPFoutputStruct.P;
r_line = 800;
z_line1 = zeros(2,numel(z));
z_line2 = zeros(2,numel(z));
for i = 1:numel(z)
   z_line1(:,i) = xObsTrue(1:2,i);
   z_line2(:,i) = [cos(z(i))*r_line; sin(z(i))*r_line] + z_line1(:,i);
end

figure(); hold on;
plot(xTgtTrue(1,:), xTgtTrue(2,:),'r');
plot(xhat_tgt_MMSE(1,:), xhat_tgt_MMSE(2,:),'k');
plot(xObsTrue(1,:), xObsTrue(2,:),'b');

markerSize = 50*w(:,end)./max(w(:,end));
markerSize(markerSize==0) = eps;
alpha_ = 0.3;
scatter(xsamps_post(1,:,end)+ xObsTrue(1,end), xsamps_post(2,:,end)+xObsTrue(2,end),markerSize,'o','MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], 'MarkerFaceAlpha', alpha_, 'MarkerEdgeAlpha', alpha_);

for i = 1:numel(z)
    plot([z_line1(1,i); z_line2(1,i)], [z_line1(2,i); z_line2(2,i)],'g');
end


legend('Target Trajectory - Truth','Target Trajectory - MMSE Estimate', 'Observer Trajectory - Truth','Final Target State pdf - PF Approximation', 'Bearing Measurements','location', 'best');
title('Engagement Ground Track');
xlabel('East (m)'); ylabel('North (m)'); 
axis equal;
% ylim([0,500])


figure(); subplot(4,1,1); hold on;
plot(timeVec,xTrue(1,:)-xhat_MMSE(1,:), 'kx');
sigma_ = squeeze(sqrt(P(1,1,:)));
plot(timeVec,0+2*sigma_, 'r--');
plot(timeVec,0-2*sigma_, 'r--');
legend('MMSE Estimate Error', 'Estimate Covariance 2\sigma');
ylabel('E Position (m)');
subplot(4,1,2); hold on;
sigma_ = squeeze(sqrt(P(2,2,:)));
plot(timeVec,xTrue(2,:)-xhat_MMSE(2,:), 'kx');
plot(timeVec,0+2*sigma_, 'r--');
plot(timeVec,0-2*sigma_, 'r--');
ylabel('N Position (m)');
legend('MMSE Estimate Error', 'Estimate Covariance 2\sigma');
xlabel('Time (s)');
subplot(4,1,3); hold on;
sigma_ = squeeze(sqrt(P(3,3,:)));
plot(timeVec,xTrue(3,:)-xhat_MMSE(3,:), 'kx');
plot(timeVec,0+2*sigma_, 'r--');
plot(timeVec,0-2*sigma_, 'r--');
ylabel('N Velocity (m/s)');
legend('MMSE Estimate Error', 'Estimate Covariance 2\sigma');
xlabel('Time (s)');
subplot(4,1,4); hold on;
sigma_ = squeeze(sqrt(P(4,4,:)));
plot(timeVec,xTrue(4,:)-xhat_MMSE(4,:), 'kx');
plot(timeVec,0+2*sigma_, 'r--');
plot(timeVec,0-2*sigma_, 'r--');
ylabel('E Velocity (m/s)');
legend('MMSE Estimate Error', 'Estimate Covariance 2\sigma');
xlabel('Time (s)');
sgtitle('Particle Filter Error');
end