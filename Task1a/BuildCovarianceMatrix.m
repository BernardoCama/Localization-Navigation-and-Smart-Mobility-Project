function [R] = buildCovarianceMatrix(parameters,TYPE)
switch TYPE
    % repet(matA, x, y) repete matA x times in horiz, y times in vertical
    case 'TOA'
        R = diag(repmat(parameters.sigmaTOA.^2,1,parameters.numberOfAP));
    case 'AOA'
        R = diag(repmat(parameters.sigmaAOA.^2,1,parameters.numberOfAP));
end

end