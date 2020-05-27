function [uHat,numberOfPerformedIterations] = iterativeNLS(parameters,AP,TYPE,R,rho)
%initilization
%uHatInit = [(rand-0.5)*parameters.xmax,(rand-0.5)*parameters.ymax];
uHatInit = [0,0];
uHat = zeros(parameters.NiterMax,2);

for iter=1:parameters.NiterMax-1
    
    % uHat(0) = uHatInit = ex (33,33)
    if iter==1
        uHat(iter,:) = uHatInit;
    end
    % uHat =  uHatInit(1) uHatInit(2)
    %         0           0
    %         0           0
    %         ..          ..
    
    %% Jacobian matrix
    [ H ] = createMatrixH(parameters,uHat(iter,:),AP,TYPE);
    
    %% evaluation of observation matrix
    [ h_uhat ] = createVectorOfObservations(parameters,uHat(iter,:),AP,TYPE);
    
    %% evaluation of Delta_rho
    delta_rho = rho - h_uhat;
    
    %% evaluation of Delta_u
    %% LS
%     IHH = inv(H'*H);
%     if rcond(IHH)  <= 1e-4
%         damp_factor = 1e-4;
%         pseudoInv = inv(H'*H + damp_factor*eye(size(H'*H)))*H';
%     else
%         pseudoInv = IHH*H';
%     end
%     delta_u = pseudoInv*delta_rho';
%   delta_u = pinv(H)*delta_rho';
    
%     %% WLS
%     IHH = inv(H'*inv(Q)*H);
%     if rcond(IHH) <= 1e-6
%         damp_factor = 1e-6;
%         pseudoInv = inv(H'*inv(Q)*H + damp_factor*eye(size(H'*inv(Q)*H)))*H';
%     else
%         pseudoInv = IHH*H';
%     end
    delta_u = inv(H' * inv(R) * H) * H' * inv(R) * delta_rho';
    
    %% update the estimate
    
    % + 0.5*delta_u' 
    uHat(iter+1,:) = uHat(iter,:) + 0.5*delta_u';
    numberOfPerformedIterations = iter+1;
    
    %% stopping criterion
    if sum(delta_u.^2)<1e-6
        return
    end       
    
end

end