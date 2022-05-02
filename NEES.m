% Function to evaluate NEES on entire block of results
function exm = NEES(X,Xh,P)

% Size parameters
Nsim = size(X,3);
N = size(X,2)-1;

% Calculate estimation errors
err = X - Xh;

% Initialize NEES statistic
exm = NaN(1,N);
ex = NaN(1,Nsim);

% Loop through each time step
for k = 1:N
    % Loop through each simulation
    for s = 1:Nsim
        ex(s) = err(:,k+1,s)'*inv(P(:,:,k+1,s))*err(:,k+1,s);
    end
    exm(k) = mean(ex);
end