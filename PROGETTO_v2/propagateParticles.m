function [PR] = propagateParticles(parameters,PR, v_, sigmaPropagation)

%sigmaPropagation = 3; %m

%motion model M1
%PR.samples = PR.samples + parameters.samplingTime.*sigmaPropagation.*randn(2,parameters.numberOfParticles);
%PR.weights = ones(1,parameters.numberOfParticles)./parameters.numberOfParticles;


%motion model M2
PR.samples(1,:) = PR.samples(1,:) + parameters.samplingTime.*(v_(1).*ones(1,parameters.numberOfParticles)  + sigmaPropagation(1).*randn(1,parameters.numberOfParticles));
PR.samples(2,:) = PR.samples(2,:) + parameters.samplingTime.*(v_(2).*ones(1,parameters.numberOfParticles)  + sigmaPropagation(2).*randn(1,parameters.numberOfParticles));

PR.weights = ones(1,parameters.numberOfParticles)./parameters.numberOfParticles;



end

