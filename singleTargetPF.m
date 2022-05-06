function [constantVelPFoutputStruct] = singleTargetPF(constantVelPFinputStruct,target_mode)
xhat0 = constantVelPFinputStruct.xhat0;
xObsTrue = constantVelPFinputStruct.xObsTrue;
xObs0 = constantVelPFinputStruct.xObs0;
z = constantVelPFinputStruct.z;
timeVec = constantVelPFinputStruct.timeVec;
xhat_MMSE = zeros(4,numel(timeVec));
Ns = constantVelPFinputStruct.Ns;
P0 = constantVelPFinputStruct.P0;
S = chol(P0, 'lower');
xsamp0 = zeros(4,Ns);
for i_NS = 1:Ns
    v_draw = S*randn(4,1);
    xsamp0(:,i_NS) = xhat0 + v_draw;
end
w = zeros(Ns,numel(timeVec));
xsamps_post = zeros(4,Ns,numel(timeVec));
xsamps = zeros(4,Ns,numel(timeVec));
% xsamps(:,:,1) = xsamp0;
P = zeros(4,4,numel(z));
effSamples = zeros(1,numel(timeVec));


constantVelPFoutputStruct = struct();
for i_z = 1:numel(z)
    
    % handle first index
    if i_z == 1
        xsamps_current = xsamp0;
        prevXObsTrue = xObs0;
        w_prev = ones(Ns,1)./Ns;
    else
        xsamps_current = xsamps(:,:,i_z-1);
        prevXObsTrue = xObsTrue(:,i_z-1);
        w_prev = w(:,i_z-1);
    end
    
    % PF update
    [xsamps_post(:,:,i_z),w(:,i_z)] = SIRPF(w_prev,xsamps_current, z(i_z), constantVelPFinputStruct, i_z, target_mode,prevXObsTrue, xObsTrue);
    effSamples(i_z) = 1/sum(w(:,i_z).^2); % calculate current sample size
    if effSamples(i_z) < 0.9*Ns % resample if particle diversity is too low
    	[xsamps(:,:,i_z), ~, ~] = RESAMPLE(xsamps_post(:,:,i_z), w(:,i_z), constantVelPFinputStruct);
        w(:,i_z) = (1/Ns).*ones(size(w(:,i_z))); % weights are now equal
    else
        xsamps(:,:,i_z) = xsamps_post(:,:,i_z); % particle diversity is adequate, do not resample
    end
    [~,I_MAP] = max(w(:,i_z));
    xhat_MAP(:,i_z) = xsamps_post(:,I_MAP,i_z); % calculate MAP estimate
    
    % calculate MMSE estimate
    xhat_MMSE(1,i_z) = sum(xsamps_post(1,:,i_z).*w(:,i_z)');
    xhat_MMSE(2,i_z) = sum(xsamps_post(2,:,i_z).*w(:,i_z)');
    xhat_MMSE(3,i_z) = sum(xsamps_post(3,:,i_z).*w(:,i_z)');
    xhat_MMSE(4,i_z) = sum(xsamps_post(4,:,i_z).*w(:,i_z)');
    
    % calculate covariance
    P(:,:,i_z) = Ns/(Ns-1) * (w(:,i_z).*(xsamps(:,:,i_z)'-xhat_MMSE(:,i_z)'))'*(xsamps(:,:,i_z)'-xhat_MMSE(:,i_z)');
    
end

constantVelPFoutputStruct.xhat_tgt_MMSE = xhat_MMSE + xObsTrue;
constantVelPFoutputStruct.P = P;
constantVelPFoutputStruct.xhat_MMSE = xhat_MMSE;
constantVelPFoutputStruct.xsamps_post = xsamps_post;
constantVelPFoutputStruct.xsamps = xsamps;
constantVelPFoutputStruct.w = w;
constantVelPFoutputStruct.effSamples = effSamples;

function [xsamps_k,w] = SIRPF(w_prev, xsamps_prev,zk, inputStruct, I_in, target_mode,prevXObsTrue, xObsTrue)
F_t = inputStruct.F_t;      % get linear motion model STM
U = inputStruct.U;          % get observer acceleration 
Q = inputStruct.Q;          % get process noise
Gamma = inputStruct.Gamma;  % get Gammma

w = zeros(1,size(xsamps_prev,2)); % preallocate particle weights

Sw = chol(Q, 'lower'); 
w_draw = Sw*randn(2,numel(w)); % draw random accels from process noise

if strcmpi(target_mode, 'straight') 
    % generate samples using 'straight' motion model
    xsamps_k = F_t(:,:,I_in)*xsamps_prev + Gamma*w_draw - U(:,I_in); 
elseif strcmpi(target_mode, 'clockwise') || strcmpi(target_mode,'counterclockwise')
    % generate samples using 'clockwise' or 'counterclockwise' motion model
    xsamps_k = zeros(size(xsamps_prev));
    for i = 1:size(xsamps_prev,2)
        % use state estimate to calculate STM
        F_t = calcF(inputStruct, prevXObsTrue,xsamps_prev(:,i),target_mode); 
        xsamps_k(:,i) = F_t*(xsamps_prev(:,i)+prevXObsTrue) - xObsTrue(:,I_in) + Gamma*w_draw(:,i);
    end
end

w = w_prev'.*get_pygivenx(zk, xsamps_k, inputStruct); 
w = w./sum(w); % renormalize 
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

    R = inputStruct.R; % get truth measurement covariance
    err = meas - atan2(x(2,:),x(1,:)); % error is zero mean with R covariance
    sigma = sqrt(R); 
    
    local_pygivenx_out = normpdf(err,0.0,sigma); % calculate p(meas|x)
end
end

