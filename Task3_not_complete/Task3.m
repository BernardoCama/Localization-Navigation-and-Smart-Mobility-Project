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

%% load measurements Task 1_a

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
% fig = figure(); hold on
% %fig.WindowState = 'maximized';
% for a = 1:parameters.numberOfAP
% 
%     imagesc(x,y,squeeze(totalLikelihood(a,:,:))');
%     
%     plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0, 254, 207]./255,'MarkerFaceColor',[0, 254, 207]./255,'DisplayName','UE')
%     
%     colorbar;
%     
%     xlabel('[m]'), ylabel('[m]');
%     
%     xlim([parameters.xmin parameters.xmax])
%     
%     ylim([parameters.ymin parameters.ymax])
%     
%     axis equal
%     
%     title(['Likelihood with TOA and AOA',' , $\sigma $ = ',num2str(parameters.sigmaTOA),' m ',' , $\sigma $ = ',num2str(rad2deg(parameters.sigmaAOA)),' deg '],'Interpreter','Latex')
% 
%     pause
% end



%% plot ML in 3D
% fig = figure()
% %fig.WindowState = 'maximized';
% for a = 1:parameters.numberOfAP
%     
%     surf(  x, y , squeeze(totalLikelihood(a,:,:)) ),hold on
%     
%     shading flat
%     
%     colorbar;
%     
%     xlabel('[m]'), ylabel('[m]');
% 
%     title(['Likelihood with TOA and AOA 3D',' , $\sigma $ = ',num2str(parameters.sigmaTOA),' m ',' , $\sigma $ = ',num2str(rad2deg(parameters.sigmaAOA)),' deg '],'Interpreter','Latex')
%     
%     plot3( UE(:,1) , UE(:,2) ,0.01, 'o','MarkerSize',10,'MarkerEdgeColor',[0, 254, 207]./255,'MarkerFaceColor',[0, 254, 207]./255, 'DisplayName','UE')
% 
%     pause
% end


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

%% Plot AP
% fig = figure(); hold on
% 
% plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255, 'DisplayName','AP')
% 
% plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0, 254, 207]./255,'MarkerFaceColor',[0, 254, 207]./255, 'DisplayName','UE')
% 
% colorbar;
% 
% xlabel('[m]'), ylabel('[m]');
% 
% xlim([parameters.xmin parameters.xmax])
% 
% ylim([parameters.ymin parameters.ymax])
% 
% axis equal
% 
% title(['Estimated Positions of AP',' , $\sigma $ = ',num2str(parameters.sigmaTOA),' m ',' , $\sigma $ = ',num2str(rad2deg(parameters.sigmaAOA)),' deg '],'Interpreter','Latex')
% 
% pause


%% load measurements Task 1_b
TYPE = 'TOA';

load('Task1b_rhoUEAP.mat')

UE = [-500, 800 ];


%% calculate R
[ h_u ] = createVectorOfObservations(parameters,UE,APhat,TYPE);

n = rhoUEAP - h_u;

R = diag(round(var(n)));









%% load measurements Task 2

load('Task2_trajectory_GR31.mat')

%% parameters
parameters.simulationTime = 200; %s
parameters.samplingTime = 1; %s
parameters.numberTrajectory = 100;

x = zeros (parameters.numberTrajectory,  parameters.simulationTime, 4);

for i = 1:parameters.numberTrajectory
    
    x(i, : ,:) = UEtrajectory{i};
    
end

% mean over all trajectories of velocity vector
v_mean = [mean(x(:,:,3),1); mean(x(:,:,4),1)];

%% plot UE trajectories
% figure,hold on
% 
% for n=1:parameters.numberTrajectory
%     
%     hold on
%     
%     plot( x(n,:,1) , x(n,:,2))
%     
%     legend('UE')
%     
%     xlabel('[m]'), ylabel('[m]');
%     
%     title('Trajectories');
%     
%     xlim([parameters.xmin parameters.xmax])
%     
%     ylim([parameters.ymin parameters.ymax])
%     
%     axis equal
%     
%     grid on
% 
% end


%% plot statistics
% figure,hold on
% 
% plot([0:parameters.samplingTime:parameters.simulationTime-1],v_mean(1,:));
% 
% plot([0:parameters.samplingTime:parameters.simulationTime-1],v_mean(2,:));
% 
% title('Velocity');xlabel('time');ylabel('m/s')
% 
acceleration = diff(v_mean(1:2,:),1,2);
% 
% figure,hold on
% 
% plot([0:parameters.samplingTime:parameters.simulationTime-2],acceleration(1,:));
% 
% plot([0:parameters.samplingTime:parameters.simulationTime-2],acceleration(2,:));
% 
% title('Acceleration');xlabel('time');ylabel('m/s2')


% since, as we can see in the graphs, the values of velocities and
% accelerations are misleading after the 100_th sample (due to constraint
% of delimited space), we calculate the statistics on samples up to 100.

sigma_tot_velocity = std(v_mean(:,1:100),0,2);
mean_tot_velocity = mean(v_mean(:,1:100),2)

sigma_tot_acceleration = std(acceleration(:,1:100),0,2);
mean_tot_acceleration = mean(acceleration(:,1:100),2);

%{
we can clearly see that is Motion Model M2 because of:
- zero sigma_tot_velocity
- constant mean_tot_velocity (about v_x = 0.4, v_y = 0.6)
- zero sigma_tot_acceleration
- zero mean_tot_acceleration

otherwise for M3 we would have obtained:
- zero sigma_tot_velocity
- zero mean_tot_velocity
- zero sigma_tot_acceleration
- zero mean_tot_acceleration

and for M4:
- zero sigma_tot_velocity
- zero mean_tot_velocity
- zero sigma_tot_acceleration
- mean_tot_acceleration != 0
%}





%% load measurements Task 3

load('Task3_rhoUEAP_GR31.mat')




%% parameters
parameters.simulationTime = 200; %s
parameters.samplingTime = 1; %s
parameters.numberTrajectory = 100;

rho = zeros (parameters.numberTrajectory,  parameters.simulationTime, 8);

for i = 1:parameters.numberTrajectory
    
    rho(i, : ,:) = rhoUEEAP{i};
    
end




%% create prior uniform distribution
parameters.numberOfParticles = 1000;

PR = generatePriorOfParticles(parameters);

figure,hold on
plot(PR.samples(1,:),PR.samples(2,:),'.r')
plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
xlabel('[m]'), ylabel('[m]');
xlim([parameters.xmin parameters.xmax])
ylim([parameters.ymin parameters.ymax])
axis equal
grid on





%% tracking
uHat = zeros(parameters.simulationTime,4);

x_mean = [0,0,mean_tot_velocity(1), mean_tot_velocity(2)]';

% covariance of Prior
P = zeros(4,4);

% motion model M2 x_t = F * x_t-1 + L * wv_t-1
T = parameters.samplingTime;
F = [1 0 T 0; 0 1 0 T; 1 0 0 0; 0 1 0 0];
L = [T 0; 0 T; 0 0; 0 0];
Q = sqrt(sum(sigma_tot_velocity.^2)) .* (L * L'); % rivedere


fig = figure(100); hold on

% tracking or the first trajectory
for time=1:parameters.simulationTime
    
    
    % Update mean and covariance
    [H] = createMatrixH(parameters,x_mean(1:2)',APhat,TYPE);
    temp = ones(size(H, 1),1);
    H = [H temp temp] ;
    
    % Kalman Gain
    G = P * H' * inv(H * P * H' + R);
    
    x_mean = x_mean + G * ( squeeze(rho(1, time, :)) - H * x_mean);
    P = P - G * H * P;
    
    uHat(time,:) = x_mean;
    

    % Prediction
    x_mean = F * x_mean;
    P = F * P * F' + Q;
    
    
    
    hold off
    
    plot( x(1,1:time,1) , x(1,1:time,2) , '-o','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
    
    hold on
    
    plot(uHat(1:time,1),uHat(1:time,2),'-s','MarkerEdgeColor',[0, 0, 160]./255,'MarkerFaceColor',[0, 0, 160]./255)
    
    hold on
    
    plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)

    legend('UE','uHat','AP')  
    xlabel('[m]'), ylabel('[m]');
    xlim([parameters.xmin parameters.xmax])
    ylim([parameters.ymin parameters.ymax])
    axis equal
    grid on
    pause(1)





end


%% distance error
% DeltaPosition = UE(:,1:2)-uHat(:,1:2) ;
% err = sqrt( sum ( DeltaPosition.^2,2));
% figure
% plot(err)





%% COMPARISON WITH PARTICLE FILTER

%{

%% Tracking by Particle Filter
parameters.numberOfParticles = 1000;

PR = generatePriorOfParticles(parameters);

figure,hold on
plot(PR.samples(1,:),PR.samples(2,:),'.r')
plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
xlabel('[m]'), ylabel('[m]');
xlim([parameters.xmin parameters.xmax])
ylim([parameters.ymin parameters.ymax])
axis equal
grid on


%% tracking
uHat = zeros(parameters.simulationTime,2);

counter = 1;

tot_PR= zeros (2,parameters.numberOfParticles, parameters.simulationTime*2);

for time=1:parameters.simulationTime
    likelihood = ones(parameters.numberOfParticles,1);
    
    %evaluate likelihood
    for a = 1:parameters.numberOfAP 
        
        likelihood = likelihood.*evaluateLikelihoodTOA2(parameters,APhat(a,:),PR,rho(1, time, a));
    end
    %normalization
    likelihood = likelihood./sum(likelihood);
    PR.weights = likelihood';
    
    indexes = resamplingAlgorithm(PR.weights,parameters.numberOfParticles);
    PR.samples = PR.samples(1:2,indexes);
    uHat(time,:) = mean(PR.samples,2);
    
%     figure,hold on
%     plot(PR.samples(1,:),PR.samples(2,:),'.r')
%     plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
%     xlabel('[m]'), ylabel('[m]');
%     xlim([parameters.xmin parameters.xmax])
%     ylim([parameters.ymin parameters.ymax])
%     axis equal
%     grid on

    tot_PR(:,:,counter) = PR.samples;
    counter = counter + 1;

    % propagation
    PR = propagateParticles(parameters,PR);
    tot_PR(:,:,counter) = PR.samples;
    counter = counter + 1;


%     figure,hold on
%     plot(PR.samples(1,:),PR.samples(2,:),'.r')
%     plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
%     xlabel('[m]'), ylabel('[m]');
%     xlim([parameters.xmin parameters.xmax])
%     ylim([parameters.ymin parameters.ymax])
%     axis equal
%     grid on


end

%% plot UE trajectory
figure,hold on
plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
plot( x(1,:,1) , x(1,:,2) , '^b')
plot( uHat(:,1) , uHat(:,2) , '*r')
legend('AP','UE GT','UE est')
xlabel('[m]'), ylabel('[m]');
xlim([parameters.xmin parameters.xmax])
ylim([parameters.ymin parameters.ymax])
axis equal
grid on

%% distance error
DeltaPosition = x(1,:,1:2)-uHat(:,1:2) ;
err = sqrt( sum ( DeltaPosition.^2,2));
figure
plot(err)

%% plot PR evolution over time


fig = figure(100); hold on

for idx = 1:counter-1
    
    hold off
    
    if mod(idx,2)
        plot (squeeze(tot_PR(1,:,idx)), squeeze(tot_PR(2,:,idx)), '.r')
    else
        plot (squeeze(tot_PR(1,:,idx)), squeeze(tot_PR(2,:,idx)), '.b')
    end
    
    hold on
    
    
    plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
    plot(  x(1,:,1) ,  x(1,:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
    plot(uHat(1:round(idx/2),1),uHat(1:round(idx/2),2),'-s','MarkerEdgeColor',[0, 0, 160]./255,'MarkerFaceColor',[0, 0, 160]./255)
    legend('Particles','AP','UE','NLS')    
    
%   legend('AP','UE','PF')
    xlabel('[m]'), ylabel('[m]');
    xlim([parameters.xmin parameters.xmax])
    ylim([parameters.ymin parameters.ymax])
    axis equal
    grid on
    
    if mod(idx,2)
        title ('Particles evolution - update', 'Interpreter', 'Latex')
    else
        title ('Particles evolution - prediction', 'Interpreter', 'Latex')
    end   
    
    pause(1)

end
%}
    


