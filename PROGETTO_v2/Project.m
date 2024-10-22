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


% plot likelihood 2D
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

APhat = circshift(APhat,4);

%% Plot AP
fig = figure(); hold on

plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255, 'DisplayName','AP')

plot( UE(:,1) , UE(:,2) , 'o','MarkerSize',10,'MarkerEdgeColor',[0, 254, 207]./255,'MarkerFaceColor',[0, 254, 207]./255, 'DisplayName','UE')

colorbar;

xlabel('[m]'), ylabel('[m]');

xlim([parameters.xmin parameters.xmax])

ylim([parameters.ymin parameters.ymax])

axis equal

title(['Estimated Positions of AP',' , $\sigma $ = ',num2str(parameters.sigmaTOA),' m ',' , $\sigma $ = ',num2str(rad2deg(parameters.sigmaAOA)),' deg '],'Interpreter','Latex')

pause





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
figure,hold on

for n=1:parameters.numberTrajectory

    hold on

    plot( x(n,:,1) , x(n,:,2))

    legend('UE')

    xlabel('[m]'), ylabel('[m]');

    title('Trajectories');

    xlim([parameters.xmin parameters.xmax])

    ylim([parameters.ymin parameters.ymax])

    axis equal

    grid on

end



%% plot statistics of first trajectory
figure,hold on

plot([0:parameters.samplingTime:parameters.simulationTime-1],v(1,:,1));

plot([0:parameters.samplingTime:parameters.simulationTime-1],v(1,:,2));

title('Velocity');xlabel('time');ylabel('m/s')

legend('Velocity_x','Velocity_y')

acceleration = squeeze(diff(v(:,:,1:2),1,2));

figure,hold on

plot([0:parameters.samplingTime:parameters.simulationTime-2],acceleration(1,:,1));

plot([0:parameters.samplingTime:parameters.simulationTime-2],acceleration(1,:,2));

title('Acceleration');xlabel('time');ylabel('m/s2')

legend('Acceleration_x','Acceleration_y')




sigma_velocity = squeeze(std(v(1,:,:),0,2))
mean_velocity = squeeze(mean(v(1,:,:),2))

sigma_acceleration = squeeze(std(acceleration(1,:,:),0,2))
mean_acceleration = squeeze(mean(acceleration(1,:,:),2))



%% motion model statistics
sigma_tot_velocity = mean(squeeze(std(v(:,:,:),0,2)),1)

sigma_tot_acceleration = mean(squeeze(std(acceleration(:,:,:),0,2)),1)
mean_tot_acceleration = mean(squeeze(mean(acceleration(:,:,:),2)),1)




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



%% tracking M2
uHat = zeros(parameters.simulationTime,4);

% which trajectory estimate
trajectory = 1;

% mean of Prior, first measurement
parameters.NiterMax = 100;
[u_0,numberOfPerformedIterations] = iterativeNLS(parameters,APhat,TYPE,R,squeeze(rho(trajectory, 1, :))');


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

Q = mean( sigma_tot_velocity,2) .* (L * L'); 


%fig = figure(100); hold on

% tracking the first trajectory
for time=1:parameters.simulationTime
    
    
    % Update mean and covariance
    [H] = createMatrixH(parameters,x_mean(1:2)',APhat,TYPE);
    
    temp = zeros(size(H, 1),1);
    
    H = [H temp temp] ;
    
    % Kalman Gain
    G = P * H' * inv(H * P * H' + R);
    
    x_mean = x_mean + G * ( squeeze(rho(trajectory, time, :)) - sqrt(sum([x_mean(1:2)'-APhat].^2,2)));
    
    P = P - G * H * P;
    
    uHat(time,:) = x_mean;
    

    % Prediction
    
    x_mean(3) = x(trajectory,time,3);
    x_mean(4) = x(trajectory,time,4);
    
    x_mean = F * x_mean;
    
    P = F * P * F' + Q;
    
   

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
    
fixes = [uHat(:,1)-x(trajectory,:,1)', uHat(:,2)-x(trajectory,:,2)', sqrt( sum ( [uHat(:,1)-x(trajectory,:,1)', uHat(:,2)-x(trajectory,:,2)'].^2,2))];
      
sorted = sort(fixes(:,3),1,'ascend');

CEP95 = sorted(length(sorted) * 0.95)

axisTheta = 0:.01:2*pi;
xunit = CEP95 * cos(axisTheta);
yunit = CEP95 * sin(axisTheta);
plot(xunit, yunit);

[in,on] = inpolygon(fixes(:,1),fixes(:,2),xunit,yunit);

plot(fixes(in,1),fixes(in,2),'r+' ,'DisplayName','IN Points'); % points inside
plot(fixes(~in,1),fixes(~in,2),'bo','DisplayName','OUT Points'); % points outside

legend('Ellipse','Fixes IN','Fixes OUT') 

title (['M2 CEP95 = ',num2str(CEP95)])

xlabel('[m]'), ylabel('[m]');

axis equal

grid on


figure,hold on

plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)

plot( x(trajectory,:,1) , x(trajectory,:,2) , '-^b')

plot( uHat(:,1) , uHat(:,2) , '-*r')

legend('AP','UE','uHat') 

title ('Task 3 Result M2')

xlabel('[m]'), ylabel('[m]');

xlim([parameters.xmin parameters.xmax])

ylim([parameters.ymin parameters.ymax])

axis equal

grid on

% pause()


%% distance error M2
DeltaPosition = squeeze(x(trajectory,:,1:2))-uHat(:,1:2) ;

err = sqrt( sum ( DeltaPosition.^2,2));

figure

plot(err)

ylim([0 8])

title ('Distance error M2')

mean_error = mean(err)




%% tracking M3
uHat = zeros(parameters.simulationTime ,4);

% which trajectory estimate
trajectory = 1;

% mean of Prior, first measurement
parameters.NiterMax = 100;
[u_0,numberOfPerformedIterations] = iterativeNLS(parameters,APhat,TYPE,R,squeeze(rho(trajectory, 1, :))');


% if we knew the velocity
 x_mean = [u_0(numberOfPerformedIterations,1),u_0(numberOfPerformedIterations, 2),0, 0]';


% covariance of Prior
P = zeros(4,4);

% motion model M3 x_t = F * x_t-1 + L * wa_t-1
% ux_t = ux_t-1 * 1 + ux_t-1 * 0 + vx_t-1 * T + vy_t-1 * 0   +  wvx_t-1 * T^2/2 + wvx_t-1 * 0
% uy_t = ux_t-1 * 0 + ux_t-1 * 1 + vx_t-1 * 0 + vy_t-1 * T   +  wvx_t-1 * 0 + wvx_t-1 * T^2/2
% vx_t = ux_t-1 * 0 + ux_t-1 * 0 + vx_t-1 * 1 + vy_t-1 * 0   +  wvx_t-1 * T + wvx_t-1 * 0
% vy_t = ux_t-1 * 0 + ux_t-1 * 0 + vx_t-1 * 0 + vy_t-1 * 1   +  wvx_t-1 * 0 + wvx_t-1 * T

T = parameters.samplingTime;

F = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];

L = [T^2/2 0; 0 T^2/2; T 0; 0 T];

Q = mean( sigma_tot_acceleration ,2) .* (L * L');


%fig = figure(100); hold on

% tracking the first trajectory
for time=1:parameters.simulationTime 
    
    
    % Update mean and covariance
    [H] = createMatrixH(parameters,x_mean(1:2)', APhat ,TYPE);
    
    temp = zeros(size(H, 1),1);
    
    H = [H temp temp] ;
    
    % Kalman Gain
    G = P * H' * inv(H * P * H' + R);
    
    x_mean = x_mean + G * ( squeeze(rho(trajectory, time, :)) - sqrt(sum([x_mean(1:2)'-APhat].^2,2)) );
    
    P = P - G * H * P;
    
    uHat(time,:) = x_mean;
    

    % Prediction
    x_mean = F * x_mean;
    
    P = F * P * F' + Q;
    
    
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
    
fixes = [uHat(:,1)-x(trajectory,:,1)', uHat(:,2)-x(trajectory,:,2)', sqrt( sum ( [uHat(:,1)-x(trajectory,:,1)', uHat(:,2)-x(trajectory,:,2)'].^2,2))];
      
sorted = sort(fixes(:,3),1,'ascend');

CEP95 = sorted(length(sorted) * 0.95)

axisTheta = 0:.01:2*pi;
xunit = CEP95 * cos(axisTheta);
yunit = CEP95 * sin(axisTheta);
plot(xunit, yunit);

[in,on] = inpolygon(fixes(:,1),fixes(:,2),xunit,yunit);

plot(fixes(in,1),fixes(in,2),'r+' ,'DisplayName','IN Points'); % points inside
plot(fixes(~in,1),fixes(~in,2),'bo','DisplayName','OUT Points'); % points outside

legend('Ellipse','Fixes IN','Fixes OUT') 

title (['M3 CEP95 = ',num2str(CEP95)])

xlabel('[m]'), ylabel('[m]');

axis equal

grid on


figure,hold on

plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)

plot( x(trajectory,:,1) , x(trajectory,:,2) , '-^b')

plot( uHat(:,1) , uHat(:,2) , '-*r')

legend('AP','UE','uHat') 

title ('Task 3 Result M3')

xlabel('[m]'), ylabel('[m]');

xlim([parameters.xmin parameters.xmax])

ylim([parameters.ymin parameters.ymax])

axis equal

grid on

% pause()


%% distance error M3
DeltaPosition = squeeze(x(trajectory,:,1:2))-uHat(:,1:2) ;

err = sqrt( sum ( DeltaPosition.^2,2));

figure

plot(err)

ylim([0 8])

title ('Distance error M3')

mean_error = mean(err)




%% COMPARISON WITH PARTICLE FILTER

%% Tracking by Particle Filter
parameters.numberOfParticles = 1000;


parameters.NiterMax = 100;
[u_0,numberOfPerformedIterations] = iterativeNLS(parameters,APhat,TYPE,R,squeeze(rho(trajectory, 1, :))');

u_0 = [u_0(numberOfPerformedIterations,1),u_0(numberOfPerformedIterations, 2)]';


PR = generatePriorOfParticles(parameters, u_0);

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

uHat(1,:) = u_0';

counter = 1;

tot_PR= zeros (2,parameters.numberOfParticles, parameters.simulationTime*2);

for time=2:parameters.simulationTime
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
    
fixes = [uHat(:,1)-x(trajectory,:,1)', uHat(:,2)-x(trajectory,:,2)', sqrt( sum ( [uHat(:,1)-x(trajectory,:,1)', uHat(:,2)-x(trajectory,:,2)'].^2,2))];
      
sorted = sort(fixes(:,3),1,'ascend');

CEP95 = sorted(length(sorted) * 0.95)

axisTheta = 0:.01:2*pi;
xunit = CEP95 * cos(axisTheta);
yunit = CEP95 * sin(axisTheta);
plot(xunit, yunit);

[in,on] = inpolygon(fixes(:,1),fixes(:,2),xunit,yunit);

plot(fixes(in,1),fixes(in,2),'r+' ,'DisplayName','IN Points'); % points inside
plot(fixes(~in,1),fixes(~in,2),'bo','DisplayName','OUT Points'); % points outside

legend('Ellipse','Fixes IN','Fixes OUT') 

title (['M2 CEP95 = ',num2str(CEP95)])

xlabel('[m]'), ylabel('[m]');

axis equal

grid on


figure,hold on

plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)

plot( x(trajectory,:,1) , x(trajectory,:,2) , '-^b')

plot( uHat(:,1) , uHat(:,2) , '-s')

legend('AP','UE GT','UE est')

title ('Task 3 Result PF')

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

title ('Distance error PF')

mean_error = mean(err)

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








%% load measurements Task 4

load('Task4_rhoUEAP_GR31.mat')




%% parameters
parameters.simulationTime = 690; %s
parameters.samplingTime = 1; %s
parameters.numberTrajectory = 1;

rho = zeros (parameters.numberTrajectory,  parameters.simulationTime, 8);

for i = 1:parameters.numberTrajectory
    
    rho(i, : ,:) = rhoUEEAP;
    
end



%% tracking M3
uHat = zeros(parameters.simulationTime ,4);

% which trajectory estimate
trajectory = 1;

% mean of Prior, first measurement
parameters.NiterMax = 100;
[u_0,numberOfPerformedIterations] = iterativeNLS(parameters,APhat,TYPE,R,squeeze(rho(trajectory, 1, :))');


x_mean = [u_0(numberOfPerformedIterations,1),u_0(numberOfPerformedIterations, 2),50, 0]';

% covariance of Prior
P = zeros(4,4);

% motion model M3 x_t = F * x_t-1 + L * wa_t-1
% ux_t = ux_t-1 * 1 + ux_t-1 * 0 + vx_t-1 * T + vy_t-1 * 0   +  wvx_t-1 * T^2/2 + wvx_t-1 * 0
% uy_t = ux_t-1 * 0 + ux_t-1 * 1 + vx_t-1 * 0 + vy_t-1 * T   +  wvx_t-1 * 0 + wvx_t-1 * T^2/2
% vx_t = ux_t-1 * 0 + ux_t-1 * 0 + vx_t-1 * 1 + vy_t-1 * 0   +  wvx_t-1 * T + wvx_t-1 * 0
% vy_t = ux_t-1 * 0 + ux_t-1 * 0 + vx_t-1 * 0 + vy_t-1 * 1   +  wvx_t-1 * 0 + wvx_t-1 * T

T = parameters.samplingTime;

F = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];

L = [T^2/2 0; 0 T^2/2; T 0; 0 T];

Q = mean( sigma_tot_acceleration ,2) .* (L * L');


%fig = figure(100); hold on

% tracking the first trajectory
for time=1:parameters.simulationTime 
    
    
    % Update mean and covariance
    [H] = createMatrixH(parameters,x_mean(1:2)', APhat ,TYPE);
    
    temp = zeros(size(H, 1),1);
    
    H = [H temp temp] ;
    
    % Kalman Gain
    G = P * H' * inv(H * P * H' + R);
    
    x_mean = x_mean + G * ( squeeze(rho(trajectory, time, :)) - sqrt(sum([x_mean(1:2)'-APhat].^2,2)) );
    
    P = P - G * H * P;
    
    uHat(time,:) = x_mean;
    

    % Prediction
    x_mean = F * x_mean;
    
    P = F * P * F' + Q;
    
    
%     hold off    
%     
%     plot(uHat(1:time,1),uHat(1:time,2),'-s','MarkerEdgeColor',[0, 0, 160]./255,'MarkerFaceColor',[0, 0, 160]./255)
%     
%     hold on
%     
%     plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
% 
%     legend('uHat','AP')  

%     xlabel('[m]'), ylabel('[m]');

%     xlim([min(uHat(1,1),uHat(time,1))-100  max(uHat(time,1),uHat(1,1))+100])

%     ylim([min(uHat(1,2),uHat(time,2))-100  max(uHat(time,2),uHat(1,2))+100])

%     grid on

%     pause(1)


end


%% plot result
figure,hold on

plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)

plot( uHat(:,1) , uHat(:,2) , '-*r')

legend('AP','uHat') 

title ('Task 4 Result M3')

xlabel('[m]'), ylabel('[m]');

xlim([parameters.xmin parameters.xmax])

ylim([parameters.ymin parameters.ymax])

axis equal

grid on

% pause()








%% Task 5
load('Task5_rhoUEAP_GR31.mat')
% %Plot Trajectory----> NLOS
hold on
plot(rhoUEEAP(:,1),'r*-');
plot(rhoUEEAP(:,3),'b*-');
plot(rhoUEEAP(:,4),'k*-');
plot(rhoUEEAP(:,2),'y*-');
plot(rhoUEEAP(:,5),'g*-');
plot(rhoUEEAP(:,6),'m*-');
plot(rhoUEEAP(:,7),'ro-');
plot(rhoUEEAP(:,8),'bo-');


%% parameters
parameters.simulationTime = 690; %s
parameters.samplingTime = 1; %s
parameters.numberTrajectory = 1;

rho = zeros (parameters.numberTrajectory,  parameters.simulationTime, 8);

for i = 1:parameters.numberTrajectory
    
    rho(i, : ,:) = rhoUEEAP;
    
end



%% tracking M3
uHat = zeros(parameters.simulationTime ,4);

% which trajectory estimate
trajectory = 1;

% mean of Prior, first measurement
parameters.NiterMax = 100;
[u_0,numberOfPerformedIterations] = iterativeNLS(parameters,APhat,TYPE,R,squeeze(rho(trajectory, 1, :))');


x_mean = [u_0(numberOfPerformedIterations,1),u_0(numberOfPerformedIterations, 2),50, 0]';

% covariance of Prior
P = zeros(4,4);

% motion model M3 x_t = F * x_t-1 + L * wa_t-1
% ux_t = ux_t-1 * 1 + ux_t-1 * 0 + vx_t-1 * T + vy_t-1 * 0   +  wvx_t-1 * T^2/2 + wvx_t-1 * 0
% uy_t = ux_t-1 * 0 + ux_t-1 * 1 + vx_t-1 * 0 + vy_t-1 * T   +  wvx_t-1 * 0 + wvx_t-1 * T^2/2
% vx_t = ux_t-1 * 0 + ux_t-1 * 0 + vx_t-1 * 1 + vy_t-1 * 0   +  wvx_t-1 * T + wvx_t-1 * 0
% vy_t = ux_t-1 * 0 + ux_t-1 * 0 + vx_t-1 * 0 + vy_t-1 * 1   +  wvx_t-1 * 0 + wvx_t-1 * T

T = parameters.samplingTime;

F = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];

L = [T^2/2 0; 0 T^2/2; T 0; 0 T];

Q = mean( sigma_tot_acceleration ,2) .* (L * L');


%fig = figure(100); hold on

% tracking the first trajectory
for time=1:parameters.simulationTime 
    
    
    % Update mean and covariance
    [H] = createMatrixH(parameters,x_mean(1:2)', APhat ,TYPE);
    
    temp = zeros(size(H, 1),1);
    
    H = [H temp temp] ;
    
    % Kalman Gain
    G = P * H' * inv(H * P * H' + R);
    
    x_mean = x_mean + G * ( squeeze(rho(trajectory, time, :)) - sqrt(sum([x_mean(1:2)'-APhat].^2,2)) );
    
    P = P - G * H * P;
    
    uHat(time,:) = x_mean;
    

    % Prediction
    x_mean = F * x_mean;
    
    P = F * P * F' + Q;
    
    
%     hold off    
%     
%     plot(uHat(1:time,1),uHat(1:time,2),'-s','MarkerEdgeColor',[0, 0, 160]./255,'MarkerFaceColor',[0, 0, 160]./255)
%     
%     hold on
%     
%     plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
% 
%     legend('uHat','AP')  

%     xlabel('[m]'), ylabel('[m]');

%     xlim([min(uHat(1,1),uHat(time,1))-100  max(uHat(time,1),uHat(1,1))+100])

%     ylim([min(uHat(1,2),uHat(time,2))-100  max(uHat(time,2),uHat(1,2))+100])

%     grid on

%     pause(1)


end


%% plot result
figure,hold on

plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)

plot( uHat(:,1) , uHat(:,2) , '-*r')

legend('AP','uHat') 

title ('Task 5 Result M3')

xlabel('[m]'), ylabel('[m]');

xlim([parameters.xmin parameters.xmax])

ylim([parameters.ymin parameters.ymax])

axis equal

grid on

% pause()






%% load measurements Task 6

load('Task6_rhoUEAP_GR31.mat')




%% parameters
parameters.simulationTime = 690; %s
parameters.samplingTime = 1; %s
parameters.numberTrajectory = 1;

rho = zeros (parameters.numberTrajectory,  parameters.simulationTime, 8);

for i = 1:parameters.numberTrajectory
    rho(i, : ,:) = rhoUEEAP;
end

hold on
plot(rhoUEEAP(:,1),'y*-');
plot(rhoUEEAP(:,2),'m*-');
plot(rhoUEEAP(:,3),'c*-');
plot(rhoUEEAP(:,4),'r*-');
plot(rhoUEEAP(:,5),'g*-');
plot(rhoUEEAP(:,6),'b*-');
plot(rhoUEEAP(:,7),'k*-');
plot(rhoUEEAP(:,8),'m*-');
legend('yellow one track 1','magenta track 2','cyan track 3','red track 4','green track 5','blue track 6','black track 7','magenta track 8');
figure

% #1
% cycle to substitute NaN values
rho_reshape=rhoUEEAP;
irows=1;
jump=0;
for icolumn=1:8
    irows=1;
    while(irows<parameters.simulationTime)
        while((irows+jump)<parameters.simulationTime && isnan(rhoUEEAP(irows+jump,icolumn)))
            jump=jump+1;
        end
        if(jump~=0 )
            correction=linspace(rhoUEEAP(irows-1,icolumn),rhoUEEAP(irows+jump,icolumn),jump+2);%
            rho_reshape(irows:irows+jump-1,icolumn) = correction (2:jump+1);
        end
        irows=irows+jump+1;
        jump=0;
    end
end

rho_reshape(isnan(rho_reshape))=0;


% #2
rho_reshape3 = zeros (parameters.simulationTime, 8);
TF = zeros (parameters.simulationTime, 8);
for i = 1:parameters.numberOfAP
    [rho_reshape3(:,i),TF(:,i)] = fillmissing(rhoUEEAP(:,i),'spline','SamplePoints',[1:parameters.simulationTime]);
end

for ir=1:8
    for ic=1:690
           if(rho_reshape(ic,ir)==0)
              rho_reshape(ic,ir)= rho_reshape3(ic,ir);
           end
    end
end




for i = 1:parameters.numberTrajectory
    rho(i, : ,:) = rho_reshape;
end

rho2 = zeros (1,parameters.simulationTime, 8);
for i = 1:parameters.numberTrajectory
    rho2(i, : ,:) = rho_reshape3;
end

% #1
hold on
plot(rho_reshape(:,1),'b*-');
plot(rhoUEEAP(:,1),'r*-');
legend('blue default 1','red fill linear 2');
title('Linear interpolation')

% #2
figure
hold on
plot(rho_reshape3(:,1),'m*-');
plot(rhoUEEAP(:,1),'g*-');
legend('blue default 1','magenta fill spline 2');
title('Cubic interpolation')


%% tracking M3 LINEAR INTERPOLATION
uHat = zeros(parameters.simulationTime ,4);

% which trajectory estimate
trajectory = 1;

% mean of Prior, first measurement
parameters.NiterMax = 100;
[u_0,numberOfPerformedIterations] = iterativeNLS_singolarity(parameters,APhat,TYPE,R,squeeze(rho(trajectory, 1, :))');


x_mean = [u_0(numberOfPerformedIterations,1),u_0(numberOfPerformedIterations, 2),50, 0]';


% covariance of Prior
P = zeros(4,4);

% motion model M3 x_t = F * x_t-1 + L * wa_t-1
% ux_t = ux_t-1 * 1 + ux_t-1 * 0 + vx_t-1 * T + vy_t-1 * 0   +  wvx_t-1 * T^2/2 + wvx_t-1 * 0
% uy_t = ux_t-1 * 0 + ux_t-1 * 1 + vx_t-1 * 0 + vy_t-1 * T   +  wvx_t-1 * 0 + wvx_t-1 * T^2/2
% vx_t = ux_t-1 * 0 + ux_t-1 * 0 + vx_t-1 * 1 + vy_t-1 * 0   +  wvx_t-1 * T + wvx_t-1 * 0
% vy_t = ux_t-1 * 0 + ux_t-1 * 0 + vx_t-1 * 0 + vy_t-1 * 1   +  wvx_t-1 * 0 + wvx_t-1 * T

T = parameters.samplingTime;

F = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];

L = [T^2/2 0; 0 T^2/2; T 0; 0 T];

Q = mean( sigma_tot_acceleration ,2) .* (L * L');


%fig = figure(100); hold on

% tracking the first trajectory
for time=1:parameters.simulationTime
    
    
    % Update mean and covariance
    [H] = createMatrixH(parameters,x_mean(1:2)', APhat ,TYPE);
    
    temp = zeros(size(H, 1),1);
    
    H = [H temp temp] ;
    [V,D] = eig(H * P * H' + R);
    if (D(1,1)<1e-4|D(2,2)<1e-4)
        diagload=1;
    else
        diagload=0;
    end
    DL=(0.0001)*(diagload*eye(8));
    
    % Kalman Gain
    G = P * H' * inv((H * P * H' + R)+diagload*DL);
    
    x_mean = x_mean + G * ( squeeze(rho(trajectory, time, :)) - sqrt(sum([x_mean(1:2)'-APhat].^2,2)) );
    
    P = P - G * H * P;
    
    uHat(time,:) = x_mean;
    
    
    % Prediction
    x_mean = F * x_mean;
    
    P = F * P * F' + Q;
    
    
    %     hold off
    %
    %     plot(uHat(1:time,1),uHat(1:time,2),'-s','MarkerEdgeColor',[0, 0, 160]./255,'MarkerFaceColor',[0, 0, 160]./255)
    %
    %     hold on
    %
    %     plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
    %
    %     legend('uHat','AP')
    
    %     xlabel('[m]'), ylabel('[m]');
    
    %     xlim([min(uHat(1,1),uHat(time,1))-100  max(uHat(time,1),uHat(1,1))+100])
    
    %     ylim([min(uHat(1,2),uHat(time,2))-100  max(uHat(time,2),uHat(1,2))+100])
    
    %     grid on
    
    %     pause(1)
    
    
end


%% plot result
figure,hold on

plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)

plot( uHat(:,1) , uHat(:,2) , '-*r')

legend('AP','uHat')

title ('Task 6 Result M3 Linear Interpolation')

xlabel('[m]'), ylabel('[m]');

xlim([parameters.xmin parameters.xmax])

ylim([parameters.ymin parameters.ymax])

axis equal

grid on

% pause()





%% tracking M3 CUBIC INTERPOLATION
uHat = zeros(parameters.simulationTime ,4);

% which trajectory estimate
trajectory = 1;

% mean of Prior, first measurement
parameters.NiterMax = 100;
[u_0,numberOfPerformedIterations] = iterativeNLS_singolarity(parameters,APhat,TYPE,R,squeeze(rho2(trajectory, 1, :))');


x_mean = [u_0(numberOfPerformedIterations,1),u_0(numberOfPerformedIterations, 2),50, 0]';


% covariance of Prior
P = zeros(4,4);

% motion model M3 x_t = F * x_t-1 + L * wa_t-1
% ux_t = ux_t-1 * 1 + ux_t-1 * 0 + vx_t-1 * T + vy_t-1 * 0   +  wvx_t-1 * T^2/2 + wvx_t-1 * 0
% uy_t = ux_t-1 * 0 + ux_t-1 * 1 + vx_t-1 * 0 + vy_t-1 * T   +  wvx_t-1 * 0 + wvx_t-1 * T^2/2
% vx_t = ux_t-1 * 0 + ux_t-1 * 0 + vx_t-1 * 1 + vy_t-1 * 0   +  wvx_t-1 * T + wvx_t-1 * 0
% vy_t = ux_t-1 * 0 + ux_t-1 * 0 + vx_t-1 * 0 + vy_t-1 * 1   +  wvx_t-1 * 0 + wvx_t-1 * T

T = parameters.samplingTime;

F = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];

L = [T^2/2 0; 0 T^2/2; T 0; 0 T];

Q = mean( sigma_tot_acceleration ,2) .* (L * L');


%fig = figure(100); hold on

% tracking the first trajectory
for time=1:parameters.simulationTime
    
    
    % Update mean and covariance
    [H] = createMatrixH(parameters,x_mean(1:2)', APhat ,TYPE);
    
    temp = zeros(size(H, 1),1);
    
    H = [H temp temp] ;
    [V,D] = eig(H * P * H' + R);
    if (D(1,1)<1e-4|D(2,2)<1e-4)
        diagload=1;
    else
        diagload=0;
    end
    DL=(0.0001)*(diagload*eye(8));
    
    % Kalman Gain
    G = P * H' * inv((H * P * H' + R)+diagload*DL);
    
    x_mean = x_mean + G * ( squeeze(rho2(trajectory, time, :)) - sqrt(sum([x_mean(1:2)'-APhat].^2,2)) );
    
    P = P - G * H * P;
    
    uHat(time,:) = x_mean;
    
    
    % Prediction
    x_mean = F * x_mean;
    
    P = F * P * F' + Q;
    
    
    %     hold off
    %
    %     plot(uHat(1:time,1),uHat(1:time,2),'-s','MarkerEdgeColor',[0, 0, 160]./255,'MarkerFaceColor',[0, 0, 160]./255)
    %
    %     hold on
    %
    %     plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
    %
    %     legend('uHat','AP')
    
    %     xlabel('[m]'), ylabel('[m]');
    
    %     xlim([min(uHat(1,1),uHat(time,1))-100  max(uHat(time,1),uHat(1,1))+100])
    
    %     ylim([min(uHat(1,2),uHat(time,2))-100  max(uHat(time,2),uHat(1,2))+100])
    
    %     grid on
    
    %     pause(1)
    
    
end


%% plot result
figure,hold on

plot( APhat(:,1) , APhat(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)

plot( uHat(:,1) , uHat(:,2) , '-*r')

legend('AP','uHat')

title ('Task 6 Result M3 Cubic Interpolation')

xlabel('[m]'), ylabel('[m]');

xlim([parameters.xmin parameters.xmax])

ylim([parameters.ymin parameters.ymax])

axis equal

grid on

% pause()




