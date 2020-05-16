function likelihood = evaluateLikelihoodAOA(parameters,rho,AP,evaluationPoint)

evaluationRho = atan2((evaluationPoint(2)-AP(2)),(evaluationPoint(1)-AP(1))); 
            
argument = rho - evaluationRho;

likelihood = 1/sqrt(2*pi*parameters.sigmaAOA.^2)*exp(-0.5* argument^2 / parameters.sigmaAOA.^2);

end