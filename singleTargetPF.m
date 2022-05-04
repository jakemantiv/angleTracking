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


constantVelPFoutputStruct = struct();
for i_z = 1:numel(z)
    
    if i_z == 1
        xsamps_current = xsamp0;
        prevXObsTrue = xObs0;
    else
        xsamps_current = xsamps(:,:,i_z-1);
        prevXObsTrue = xObsTrue(:,i_z-1);
    end
    [xsamps_post(:,:,i_z),w(:,i_z)] = SIRPF(xsamps_current, z(i_z), constantVelPFinputStruct, i_z, target_mode,prevXObsTrue, xObsTrue);
    [xsamps(:,:,i_z), ~, ~] = RESAMPLE(xsamps_post(:,:,i_z), w(:,i_z), constantVelPFinputStruct);
%     [xsamps(:,:,i_z), ~] = FAST_RESAMPLE(xsamps_post(:,:,i_z), w(:,i_z), constantVelPFinputStruct);
%     [uv, IA, IC] = unique(xsamps_post(:,:,i)','rows');
    [~,I_MAP] = max(w(:,i_z));
%     xhat_MAP(:,i) = xsamps_post(:,I_MAP,i);
    xhat_MMSE(1,i_z) = sum(xsamps_post(1,:,i_z).*w(:,i_z)');
    xhat_MMSE(2,i_z) = sum(xsamps_post(2,:,i_z).*w(:,i_z)');
    xhat_MMSE(3,i_z) = sum(xsamps_post(3,:,i_z).*w(:,i_z)');
    xhat_MMSE(4,i_z) = sum(xsamps_post(4,:,i_z).*w(:,i_z)');
    P(:,:,i_z) = cov(xsamps(:,:,i_z)');
%     cv = var(w)/(mean(w)^2);
%     effSamples(i) = Ns/(1+cv);
    
end

constantVelPFoutputStruct.xhat_tgt_MMSE = xhat_MMSE + xObsTrue;
constantVelPFoutputStruct.P = P;
constantVelPFoutputStruct.xhat_MMSE = xhat_MMSE;
constantVelPFoutputStruct.xsamps_post = xsamps_post;
constantVelPFoutputStruct.xsamps = xsamps;
constantVelPFoutputStruct.w = w;

function [xsamps_k,w] = SIRPF(xsamps_prev,zk, inputStruct, I_in, target_mode,prevXObsTrue, xObsTrue)
F_t = inputStruct.F_t;
% G = inputStruct.G;
U = inputStruct.U;
Q = inputStruct.Q;
Gamma = inputStruct.Gamma;
Sw = chol(Q, 'lower'); 

w = zeros(1,size(xsamps_prev,2));
w_draw = Sw*randn(2,numel(w));
if strcmpi(target_mode, 'straight')
    xsamps_k = F_t(:,:,I_in)*xsamps_prev + Gamma*w_draw - U(:,I_in);
elseif strcmpi(target_mode, 'clockwise') || strcmpi(target_mode,'counterclockwise')
    xsamps_k = F_t(:,:,I_in)*(xsamps_prev+prevXObsTrue) - xObsTrue(:,I_in) + Gamma*w_draw;
end
w = get_pygivenx(zk, xsamps_k, inputStruct);
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
function [xsamps_out, w_out] = FAST_RESAMPLE(xsamps_in, w_in, inputStruct)

nsampsIn = numel(w_in);
sampcdf = repmat(cumsum(w_in),[1 1]);

%draw a bunch of random #'s from U(0,1):
urands = repmat(rand(1,nsampsIn),nsampsIn,1);

%compare to CDF:
tag = (urands<sampcdf); %returns a logical array
tag = single(tag); %convert to a numeric array
tagsums = cumsum(tag,1); %tally up # of times tomcruise is less than cdf

%for each sample (column) extract the rows where the first < occurs in tag
[indsampsout,~] = find(tagsums==1);
xsamps_out = xsamps_in(:,indsampsout);
w_out = 1/nsampsIn*ones(size(w_in));
end
function [local_pygivenx_out] = get_pygivenx(meas, x, inputStruct)

    R = inputStruct.R; 
    err = meas - atan2(x(2,:),x(1,:)); % error is zero mean with R covariance
    sigma = sqrt(R);
    
    local_pygivenx_out = normpdf(err,0.0,sigma);
end
end

