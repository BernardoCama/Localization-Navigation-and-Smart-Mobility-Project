function likelihood = evaluateLikelihoodTOA(parameters,AP,PR,rho)

 
evaluationRange = sqrt(sum([PR.samples(1:2,:)'-AP].^2,2)); 
argument = rho - evaluationRange;

likelihood = 1./sqrt(2*pi*parameters.sigmaTOA.^2).*exp(-0.5.* (argument./parameters.sigmaTOA).^2);


end