clear all, clc, close all

set(0,'DefaultTextFontSize',18)
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',16)

%% scenario settings (4000x4000 m)
parameters.xmin = -2000; parameters.ymin = -2000;
parameters.xmax =  2000; parameters.ymax =  2000;

%% position of the Access Points
parameters.numberOfAP = 8;
% [ AP ] = generatePositionOfAP(parameters);

%% position of the user
UE = [0,0];

%% accuracies
parameters.sigmaTOA = 20; %m
parameters.sigmaAOA = deg2rad(10); %deg

%% create a grid of evaluation points
x = linspace(parameters.xmin,parameters.xmax,1000);
y = linspace(parameters.ymin,parameters.ymax,1000);

%% build covariance matrix
% TYPE='TOA';
% [R] = BuildCovarianceMatrix(parameters,TYPE);

%% load measurements

load('Task1a_rhoUEAP.mat')

rho=rhoUEAP;

%% compute likelihood for each AP in each evaluation point
likelihood = zeros(parameters.numberOfAP,length(x),length(y), 2);


for a = 1:parameters.numberOfAP

    % Evaluate the likelihood in each evaluation point
    for i=1:length(x)
        for j=1:length(y)

                    likelihood(a,i,j,1) = evaluateLikelihoodTOA(parameters,rho(a,1),UE,[x(i),y(j)]);

                    likelihood(a,i,j,2) = evaluateLikelihoodAOA(parameters,rho(a,2),UE,[x(i),y(j)]);

        end %j
    end %i
end %a  

%% compute "global" ML likelihood
totalLikelihood = zeros(parameters.numberOfAP, length(x),length(y));

% Multiply every likelihood(AP) and then normalize
for a = 1:parameters.numberOfAP
    
    totalLikelihood(a,:,:) = squeeze(likelihood(a,:,:,1).*likelihood(a,:,:,2));
    
    totalLikelihood(a) = totalLikelihood(a)./sum(sum(totalLikelihood(a)));
end


%% plot likelihood 2D
fig = figure(); hold on
%fig.WindowState = 'maximized';
for a = 1:parameters.numberOfAP

    imagesc(x,y,squeeze(totalLikelihood(a,:,:))');
    
    plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0, 254, 207]./255,'MarkerFaceColor',[0, 254, 207]./255,'DisplayName','UE')
    
    colorbar;
    
    xlabel('[m]'), ylabel('[m]');
    
    xlim([parameters.xmin parameters.xmax])
    
    ylim([parameters.ymin parameters.ymax])
    
    axis equal
    
    title(['Likelihood with TOA and AOA',' , $\sigma $ = ',num2str(parameters.sigmaTOA),' m ',' , $\sigma $ = ',num2str(rad2deg(parameters.sigmaAOA)),' deg '],'Interpreter','Latex')

    pause
end



%% plot ML in 3D
fig = figure()
%fig.WindowState = 'maximized';
for a = 1:parameters.numberOfAP
    
    surf(  x, y , squeeze(totalLikelihood(a,:,:)) ),hold on
    
    shading flat
    
    colorbar;
    
    xlabel('[m]'), ylabel('[m]');

    title(['Likelihood with TOA and AOA 3D',' , $\sigma $ = ',num2str(parameters.sigmaTOA),' m ',' , $\sigma $ = ',num2str(rad2deg(parameters.sigmaAOA)),' deg '],'Interpreter','Latex')
    
    plot3( UE(:,1) , UE(:,2) ,0.01, 'o','MarkerSize',10,'MarkerEdgeColor',[0, 254, 207]./255,'MarkerFaceColor',[0, 254, 207]./255, 'DisplayName','UE')

    pause
end


%% evaluate the ML
maxValue = zeros (parameters.numberOfAP,1);

APhat = zeros (parameters.numberOfAP,2);

for a = 1:parameters.numberOfAP
    
    temp = squeeze(totalLikelihood(a,:,:));
    
    maxValue(a) = max(temp(:));
    
    [xhat yhat] = find(temp == maxValue(a));
    
    APhat(a,1) = round(x(xhat));
    
    APhat(a,2) = round(y(yhat));

end


fig = figure(); hold on

plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255, 'DisplayName','AP')

plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0, 254, 207]./255,'MarkerFaceColor',[0, 254, 207]./255, 'DisplayName','UE')

colorbar;

xlabel('[m]'), ylabel('[m]');

xlim([parameters.xmin parameters.xmax])

ylim([parameters.ymin parameters.ymax])

axis equal

title(['Estimated Positions of AP',' , $\sigma $ = ',num2str(parameters.sigmaTOA),' m ',' , $\sigma $ = ',num2str(rad2deg(parameters.sigmaAOA)),' deg '],'Interpreter','Latex')




