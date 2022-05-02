clear; close all;
saveFigs = false;
rng(12);

% favorable geometry
obsPos0 = [-300; 10]; 
obsVel0 = [20; 1];

% unfavorable geometry
% obsPos0 = [-10; 10]; 
% obsVel0 = [1; 0];


tgtPos0 = [50; 90];
tgtVel0 = [1; 6];

xTgt0 = [tgtPos0; tgtVel0];
xObs0 = [obsPos0; obsVel0];

xTrue0 = xTgt0 - xObs0;

dt = 5; % sampling interval seconds
maxTime = 50;
Gamma = [dt^2/2, 0;
         0,  dt^2/2;
         dt, 0
         0, dt];
F = [1, 0, dt, 0;
     0, 1, 0, dt;
     0, 0, 1, 0; 
     0, 0, 0, 1];
 
timeVec = dt:dt:maxTime;

Q = 0.001.*eye(2);

R = deg2rad(0.01);

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



xhat_MMSE = zeros(size(xTrue));
inputStruct.R = R;
inputStruct.Q = Q; 
inputStruct.F = F;
% inputStruct.U = U;
inputStruct.Gamma = Gamma; 
Ns = 1000;

P0 = [100, 0, 0, 0; 
      0, 100, 0, 0;
      0, 0, 10, 0;
      0, 0, 0, 10];
S = chol(P0, 'lower');
xsamp = zeros(4,Ns);
for i = 1:Ns
    v_draw = S*randn(4,1);
    xsamp(:,i) = xTrue0 + v_draw;
end
w = zeros(Ns,numel(timeVec));
xsamps_post = zeros(4,Ns,numel(timeVec));
P = zeros(4,4,numel(z));
for i = 1:numel(z)

    [xsamps_post(:,:,i),w(:,i)] = SIRPF(xsamp, z(i), inputStruct, i);
    [xsamp, ~, ~] = RESAMPLE(xsamps_post(:,:,i), w(:,i), inputStruct);
%     [uv, IA, IC] = unique(xsamps_post(:,:,i)','rows');
    [~,I_MAP] = max(w(:,i));
%     xhat_MAP(:,i) = xsamps_post(:,I_MAP,i);
    xhat_MMSE(1,i) = sum(xsamps_post(1,:,i).*w(:,i)');
    xhat_MMSE(2,i) = sum(xsamps_post(2,:,i).*w(:,i)');
    xhat_MMSE(3,i) = sum(xsamps_post(3,:,i).*w(:,i)');
    xhat_MMSE(4,i) = sum(xsamps_post(4,:,i).*w(:,i)');
    P(:,:,i) = cov(xsamp');
%     cv = var(w)/(mean(w)^2);
%     effSamples(i) = Ns/(1+cv);
    
end

xhat_tgt_MMSE = xhat_MMSE + xObsTrue;

r_line = 800;
z_line1 = zeros(2,numel(z));
z_line2 = zeros(2,numel(z));
for i = 1:numel(z)
   z_line1(:,i) = xObsTrue(1:2,i);
   z_line2(:,i) = [cos(z(i))*r_line; sin(z(i))*r_line] + z_line1(:,i);
end


figure(); hold on;
plot3(xTgtTrue(1,:), xTgtTrue(2,:), timeVec,'r');
plot3(xhat_tgt_MMSE(1,:), xhat_tgt_MMSE(2,:), timeVec,'k');

plot3(xObsTrue(1,:), xObsTrue(2,:), timeVec,'b');
for i = 1:numel(z)
    plot3([z_line1(1,i); z_line2(1,i)], [z_line1(2,i); z_line2(2,i)], [timeVec(i); timeVec(i)],'g');
end

markerSize = 50*w(:,end)./max(w(:,end));
alpha_ = 0.5;
scatter3(xsamps_post(1,:,end)+ xObsTrue(1,end), xsamps_post(2,:,end)+xObsTrue(2,end),timeVec(end).*ones(Ns,1),markerSize,'o','MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], 'MarkerFaceAlpha', alpha_, 'MarkerEdgeAlpha', alpha_);
legend('Target Trajectory - Truth','Target Trajectory - Estimate', 'Observer Trajectory - Truth', 'Bearing Measurements');
title('Engagement Ground Track');
xlabel('East (m)'); ylabel('North (m)'); 
axis equal; 
% ylim([0,500])


figure(); hold on;
plot(xTgtTrue(1,:), xTgtTrue(2,:),'r');
plot(xhat_tgt_MMSE(1,:), xhat_tgt_MMSE(2,:),'k');
plot(xObsTrue(1,:), xObsTrue(2,:),'b');

markerSize = 50*w(:,end)./max(w(:,end));
alpha_ = 0.3;
scatter(xsamps_post(1,:,end)+ xObsTrue(1,end), xsamps_post(2,:,end)+xObsTrue(2,end),markerSize,'o','MarkerFaceColor',[0.8500, 0.3250, 0.0980],'MarkerEdgeColor',[0.8500, 0.3250, 0.0980], 'MarkerFaceAlpha', alpha_, 'MarkerEdgeAlpha', alpha_);

for i = 1:numel(z)
    plot([z_line1(1,i); z_line2(1,i)], [z_line1(2,i); z_line2(2,i)],'g');
end


legend('Target Trajectory - Truth','Target Trajectory - MMSE Estimate', 'Observer Trajectory - Truth','Final Target State pdf - PF Approximation', 'Bearing Measurements','location', 'best');
title('Engagement Ground Track');
xlabel('East (m)'); ylabel('North (m)'); 
axis equal;
ylim([0,500])


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


function [obsStateNext] = constantVelObsModel(obsStatePrev,dt)
    F_obs = [1, 0, dt, 0;
         0, 1, 0, dt;
         0, 0, 1, 0; 
         0, 0, 0, 1];
    obsStateNext = F_obs*obsStatePrev;
end
function [xsamps_k,w] = SIRPF(xsamps_prev,zk, inputStruct, I_in)
F = inputStruct.F;
% G = inputStruct.G;
% U = inputStruct.U;
Q = inputStruct.Q;
R = inputStruct.R;
Gamma = inputStruct.Gamma;
Sw = chol(Q, 'lower'); 

% pygivenx = zeros(Nx,1);
xsamps_k = zeros(size(xsamps_prev));
w = zeros(1,size(xsamps_prev,2));
for i = 1:size(xsamps_prev,2)
%    xsamp_vec = zeros(Nx,1);
%    xsamp_vec(xsamps_prev(i)) = 1;
%    xsamps_k(i) = randsample(Nx, 1, true,T*xsamp_vec);
   w_draw = Sw*randn(2,1);
%    xsamps_k(:,i) = F*xsamps_prev(:,i) + G*U(I_in) + w_draw;
   xsamps_k(:,i) = F*xsamps_prev(:,i) + Gamma*w_draw;

   %    xsamps_k(:,i) = F*x(:,i-1) + G*U(i) + w_draw;
   w(i) = get_pygivenx(zk, xsamps_k(:,i), inputStruct);
end
w = w./sum(w);




end

function [xsamps_out, w_out, parent_out] = RESAMPLE(xsamps_in, w_in, inputStruct)
w_out = zeros(size(w_in));
xsamps_out = zeros(size(xsamps_in));
parent_out = zeros(size(w_in));

i = 1;
c = cumsum(w_in);
N = numel(w_in);
u1 = 0 + (1/N)*rand(1); % u1 ~ U[0, 1/N]
for j = 1:N
   uj = u1 + (j-1)/N;
   while uj > c(i)
       i = i+1;
   end
   xsamps_out(:,j) = xsamps_in(:,i);
   w_out(j) = 1/N;
   parent_out(j) = i; 
end
end

function [local_pygivenx_out] = get_pygivenx(meas, x, inputStruct)
    R = inputStruct.R; 
    mean = atan2(x(2),x(1));
    sigma = sqrt(R);
    
    local_pygivenx_out = normpdf(meas,mean,sigma);
    
end
% function [tgtStateNext] = constantVelTgtModel(tgtStatePrev,dt)
%     F_tgt = [1, 0, dt, 0;
%          0, 1, 0, dt;
%          0, 0, 1, 0; 
%          0, 0, 0, 1];
%     tgtStateNext = F_tgt*tgtStatePrev;
% end