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

APhat = circshift(APhat,4);

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

UE = [500, -800 ];


%% calculate R
[ h_u ] = createVectorOfObservations(parameters,UE,APhat,TYPE);

n = rhoUEAP - h_u;

R = diag(round(var(n)))









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

% extract velocity
v = x(:,:,3:4);

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



%% plot statistics of first trajectory
figure,hold on

plot([0:parameters.samplingTime:parameters.simulationTime-1],v(1,:,1));

plot([0:parameters.samplingTime:parameters.simulationTime-1],v(1,:,2));

title('Velocity');xlabel('time');ylabel('m/s')

acceleration = squeeze(diff(v(:,:,1:2),1,2));

figure,hold on

plot([0:parameters.samplingTime:parameters.simulationTime-2],acceleration(1,:,2));

plot([0:parameters.samplingTime:parameters.simulationTime-2],acceleration(1,:,2));

title('Acceleration');xlabel('time');ylabel('m/s2')



sigma_velocity = squeeze(std(v(1,:,:),0,2))
mean_velocity = squeeze(mean(v(1,:,:),2))

sigma_acceleration = squeeze(std(acceleration(1,:,:),0,2))
mean_acceleration = squeeze(mean(acceleration(1,:,:),2))



%% motion model statistics
sigma_tot_velocity = mean(squeeze(std(v(:,:,:),0,2)),1)

sigma_tot_acceleration = mean(squeeze(std(acceleration(:,:,:),0,2)),1)
mean_tot_acceleration = mean(squeeze(mean(acceleration(:,:,:),2)),1)

%{
RIVEDERE SBAGLIATO
we can clearly see that is Motion Model M2 because of:
- zero sigma_tot_velocity
- constant mean_tot_velocity (about v_x = 1.387, v_y = 5.635)
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



%% tracking
uHat = zeros(parameters.simulationTime,4);

% which trajectory estimate
trajectory = 1;

% mean of Prior, first measurement
parameters.NiterMax = 100;
[u_0,numberOfPerformedIterations] = iterativeNLS(parameters,APhat,TYPE,R,squeeze(rho(trajectory, 1, :))');

% if we don't knew the velocity
% x_mean = [u_0(numberOfPerformedIterations,1),u_0(numberOfPerformedIterations, 2),mean(x(trajectory,:,3),2), mean(x(trajectory,:,4),2)]';

% if we knew the velocity
x_mean = [u_0(numberOfPerformedIterations,1),u_0(numberOfPerformedIterations, 2),x(trajectory,1,3), x(trajectory,1,4)]';




% covariance of Prior
P = zeros(4,4);

% motion model M2 x_t = F * x_t-1 + L * wv_t-1
% ux_t = ux_t-1 * 1 + ux_t-1 * 0 + vx_t-1 * T + vy_t-1 * 0   +  wvx_t-1 * T + wvx_t-1 * 0
% uy_t = ux_t-1 * 0 + ux_t-1 * 1 + vx_t-1 * 0 + vy_t-1 * T   +  wvx_t-1 * 0 + wvx_t-1 * T
% vx_t = ux_t-1 * 0 + ux_t-1 * 0 + vx_t-1 * 1 + vy_t-1 * 0   +  wvx_t-1 * 0 + wvx_t-1 * 0
% vy_t = ux_t-1 * 0 + ux_t-1 * 0 + vx_t-1 * 0 + vy_t-1 * 1   +  wvx_t-1 * 0 + wvx_t-1 * 0

T = parameters.samplingTime;

F = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];

L = [T 0; 0 T; 0 0; 0 0];

Q = 0 .* (L * L'); % (CAMBIARE DA TASK2)


%fig = figure(100); hold on

% tracking the first trajectory
for time=1:parameters.simulationTime
    
    
    % Update mean and covariance
    [H] = createMatrixH(parameters,x_mean(1:2)',APhat,TYPE);
    
    temp = ones(size(H, 1),1);
    
    H = [H temp temp] ;
    
    % Kalman Gain
    G = P * H' * inv(H * P * H' + R);
    
    x_mean = x_mean + G * ( squeeze(rho(trajectory, time, :)) - H * x_mean);
    
    P = P - G * H * P;
    
    uHat(time,:) = x_mean;
    

    % Prediction
    
    % if we knew the velocity
    x_mean(3) = x(trajectory,time,3);
    x_mean(4) = x(trajectory,time,4);
    
    x_mean = F * x_mean;
    
    P = F * P * F' + Q;
    
    
    % P aumenta sempre di piu quando ci si somma Q, 
    % non diminusce molto quando ci si sottrae - G * H * P;
    
    % x_mean = F * x_mean; aumenta del giusto
    
    
%     hold off
%     
%     plot( x(trajectory,1:time,1) , x(trajectory,1:time,2) , '-o','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
%     
%     hold on
%     
%     plot(uHat(1:time,1),uHat(1:time,2),'-s','MarkerEdgeColor',[0, 0, 160]./255,'MarkerFaceColor',[0, 0, 160]./255)
%     
%     hold on
%     
%     plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
% 
%     legend('UE','uHat','AP')  

%     xlabel('[m]'), ylabel('[m]');

%     xlim([min(x(trajectory,1,1),uHat(1,1))-100  max(x(trajectory,time,1),uHat(time,1))+100])

%     ylim([min(x(trajectory,1,2),uHat(1,2))-100  max(x(trajectory,time,2),uHat(time,2))+100])

%     grid on

%     pause(1)


end


%% plot result
figure,hold on

plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)

plot( x(trajectory,:,1) , x(trajectory,:,2) , '-^b')

plot( uHat(:,1) , uHat(:,2) , '-*r')

legend('AP','UE','uHat') 

xlabel('[m]'), ylabel('[m]');

xlim([parameters.xmin parameters.xmax])

ylim([parameters.ymin parameters.ymax])

axis equal

grid on

% pause()


%% distance error
DeltaPosition = squeeze(x(1,:,1:2))-uHat(:,1:2) ;

err = sqrt( sum ( DeltaPosition.^2,2));

figure

plot(err)

title ('Distance error')




%% COMPARISON WITH PARTICLE FILTER

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
        
        likelihood = likelihood.*evaluateLikelihoodTOA2(parameters,APhat(a,:),PR,rho(trajectory, time, a));
    
    end
    
    %normalization
    likelihood = likelihood./sum(likelihood);
    
    PR.weights = likelihood';
    
    indexes = resamplingAlgorithm(PR.weights,parameters.numberOfParticles);
    
    PR.samples = PR.samples(1:2,indexes);
    
    uHat(time,:) = mean(PR.samples,2);
    
%     figure,hold on

%     plot(PR.samples(1,:),PR.samples(2,:),'.r')

%     plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)

%     xlabel('[m]'), ylabel('[m]');

%     xlim([parameters.xmin parameters.xmax])

%     ylim([parameters.ymin parameters.ymax])

%     axis equal

%     grid on

    tot_PR(:,:,counter) = PR.samples;
    
    counter = counter + 1;

    % propagation
    PR = propagateParticles(parameters,PR, [mean(x(trajectory,:,3),2) mean(x(trajectory,:,4),2)], sigma_tot_velocity);
    
    tot_PR(:,:,counter) = PR.samples;
    
    counter = counter + 1;


%     figure,hold on

%     plot(PR.samples(1,:),PR.samples(2,:),'.r')

%     plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)

%     xlabel('[m]'), ylabel('[m]');

%     xlim([parameters.xmin parameters.xmax])

%     ylim([parameters.ymin parameters.ymax])

%     axis equal

%     grid on


end

%% plot UE trajectory
figure,hold on

plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)

plot( x(trajectory,:,1) , x(trajectory,:,2) , '-^b')

plot( uHat(:,1) , uHat(:,2) , '-s')

legend('AP','UE GT','UE est')

xlabel('[m]'), ylabel('[m]');

xlim([parameters.xmin parameters.xmax])

ylim([parameters.ymin parameters.ymax])

axis equal

grid on


%% distance error
DeltaPosition = squeeze(x(trajectory,:,1:2))-uHat(:,1:2) ;

err = sqrt( sum ( DeltaPosition.^2,2));

figure

plot(err)

title ('Distance error')

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
    
    plot(  x(trajectory,1:round(idx/2),1) ,  x(trajectory,1:round(idx/2),2) , 'o','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
    
    plot(uHat(1:round(idx/2),1),uHat(1:round(idx/2),2),'-s','MarkerEdgeColor',[0, 0, 160]./255,'MarkerFaceColor',[0, 0, 160]./255)
    
    legend('Particles','AP','UE','NLS')    
    
    xlabel('[m]'), ylabel('[m]');
    
    xlim([min(x(trajectory,1,1),uHat(1,1))-100  max(x(trajectory,round(idx/2),1),uHat(round(idx/2),1))+100])
    
    ylim([min(x(trajectory,1,2),uHat(1,2))-100  max(x(trajectory,round(idx/2),2),uHat(round(idx/2),2))+100])
    
    %axis equal
    
    grid on
    
    if mod(idx,2)
        title ('Particles evolution - update', 'Interpreter', 'Latex')
    else
        title ('Particles evolution - prediction', 'Interpreter', 'Latex')
    end   
    
    pause(1)

end





