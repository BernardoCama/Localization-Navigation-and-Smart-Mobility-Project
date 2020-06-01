function [PR] = propagateParticles(parameters,PR, v_, sigmaPropagation)

sigmaPropagation = 3; %m

%motion model
PR.samples = PR.samples + parameters.samplingTime.*sigmaPropagation.*randn(2,parameters.numberOfParticles);
% PR.weights = ones(1,parameters.numberOfParticles)./parameters.numberOfParticles;

%v_ = [mean(x(1,:,3),2) mean(x(1,:,4),2)]
%sigmaPropagation = sigma_tot_velocity;

%motion model M2
%PR.samples = PR.samples + parameters.samplingTime.*[v_(1).*ones(1,parameters.numberOfParticles) ;v_(2).*ones(1,parameters.numberOfParticles)] +parameters.samplingTime.*[sigmaPropagation(1).*randn(1,parameters.numberOfParticles) ; sigmaPropagation(2).*randn(1,parameters.numberOfParticles)];

%PR.samples = PR.samples + v + parameters.samplingTime.*[sigmaPropagation(1).*randn(1,parameters.numberOfParticles) ; sigmaPropagation(2).*randn(1,parameters.numberOfParticles)];

PR.weights = ones(1,parameters.numberOfParticles)./parameters.numberOfParticles;

end

