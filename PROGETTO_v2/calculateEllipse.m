function  [ellipsePoints] = calculateEllipse ( parameters, H, R , UE , AP , TYPE , k) 
% input k is the k-sigma ellipse value (e.g., k=3 sigma)

% singular value decomposition
C = inv(H'*inv(R)*H);

%rcond(C) returns an estimate for the reciprocal condition of C in 1-norm.
%If C is well conditioned, rcond(C) is near 1.0.
%If C is badly conditioned, rcond(C) is near 0.
if rcond(C) <= 1e-4
    damp_factor = 1e-4;
    C = inv(H'*inv(R)*H + damp_factor*eye(size(H'*inv(R)*H)));
end
[U, S, ~] = svd(C); %[U,S,V] = svd(C) -> C = U*S*V' 
                    % ~ -> don't care
                    
% Calculating the SVD consists of finding the eigenvalues and eigenvectors of C*C' and C*C'.
% The eigenvectors of C*C' make up the columns of V , the eigenvectors of C*C'  make up the columns of U. 
% Also, the singular values in S are square roots of eigenvalues from C*C' or C*C'.
%
%
S = sqrt(S);  % eigenvector -> lambda1, lambda2 = sigma_H1, sigma_H2      
[~,pos_max] = max(diag(S));

a = k*max(max(S));  % Major semiaxis -> k*sigma_H1 
b = k*min(max(S));  % Minor semiaxis -> k*sigma_H2
theta = atan2(U(2,pos_max),U(1,pos_max));  % orientation of the ellipse
                                           % U are the eigenvector

% ellipse equation
axisTheta = 0:.01:2*pi;
ellipsePoints(1,:) = a*cos(axisTheta)*cos(theta) - b*sin(axisTheta)*sin(theta) + UE(1);
ellipsePoints(2,:) = a*cos(axisTheta)*sin(theta) + b*sin(axisTheta)*cos(theta) + UE(2);

% Plot ellipse
% plotellipse(parameters,AP,ellipsePoints,UE,TYPE)


end