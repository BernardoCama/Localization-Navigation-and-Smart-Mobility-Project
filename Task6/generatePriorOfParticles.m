function PR = generatePriorOfParticles(parameters)

NP = parameters.numberOfParticles;

PR.samples = zeros (2, NP);
PR.weights = zeros (1, NP);

% prior pdf
PR.samples(1,:) = (rand([1, NP]) - 0.5).*2.*parameters.xmax;
PR.samples(2,:) = (rand([1, NP]) - 0.5).*2.*parameters.ymax;
PR.weights = 1./NP * ones(1,NP);

end

