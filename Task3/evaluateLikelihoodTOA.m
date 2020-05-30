function likelihood = evaluateLikelihoodTOA(parameters,rho,AP,evaluationPoint)

evaluationRho = sqrt(sum([evaluationPoint-AP].^2,2)); % h(u)
argument = rho - evaluationRho ;

likelihood = 1/sqrt(2*pi*parameters.sigmaTOA.^2)*exp(-0.5*(argument/parameters.sigmaTOA)^2);


end