function [PR] = propagateParticles(parameters,PR)

sigmaPropagation = 3; %m

%motion model M1
PR.samples = PR.samples + parameters.samplingTime.*sigmaPropagation.*randn(2,parameters.numberOfParticles);
PR.weights = ones(1,parameters.numberOfParticles)./parameters.numberOfParticles;
end