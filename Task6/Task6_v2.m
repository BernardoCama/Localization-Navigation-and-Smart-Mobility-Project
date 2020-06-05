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

plot([0:parameters.samplingTime:parameters.simulationTime-2],acceleration(1,:,1));

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
            rho_reshape(irows:irows+jump-1,icolumn)=linspace(rhoUEEAP(irows-1,icolumn),rhoUEEAP(irows+jump,icolumn),jump);%
        end
        irows=irows+jump+1;
        jump=0;
    end
end

rho_reshape(isnan(rho_reshape))=0;
for i = 1:parameters.numberTrajectory
    rho(i, : ,:) = rho_reshape;
end


% #2
rho_reshape2 = zeros (parameters.simulationTime, 8);
TF = zeros (parameters.simulationTime, 8);
for i = 1:parameters.numberOfAP
    [rho_reshape2(:,i),TF(:,i)] = fillmissing(rhoUEEAP(:,i),'linear','SamplePoints',[1:parameters.simulationTime]);
end


% #3
rho_reshape3 = zeros (parameters.simulationTime, 8);
TF = zeros (parameters.simulationTime, 8);
for i = 1:parameters.numberOfAP
    [rho_reshape3(:,i),TF(:,i)] = fillmissing(rhoUEEAP(:,i),'spline','SamplePoints',[1:parameters.simulationTime]);
end


% #1
hold on
plot(rho_reshape(:,1),'b*-');
plot(rhoUEEAP(:,1),'r*-');
legend('blue default 1','red algorithm 2');

% #2
figure
hold on
plot(rho_reshape2(:,1),'b*-');
plot(rhoUEEAP(:,1),'g*-');
legend('blue default 1','green fill linear 2');

% #3
figure
hold on
plot(rho_reshape3(:,1),'m*-');
plot(rhoUEEAP(:,1),'g*-');
legend('blue default 1','magenta fill spline 2');


% % plot(rhocorrect(:,3),'b*-');
% % plot(rhocorrect(:,4),'k*-');
% % plot(rhocorrect(:,2),'y*-');
% % plot(rhocorrect(:,5),'g*-');
% % plot(rhocorrect(:,6),'m*-');
% % plot(rhocorrect(:,7),'m*-');
% % plot(rhocorrect(:,8),'m*-');
%legend('red one track 1','yellow track 2','blue track 3','black track 4','green track 5','magenta track 6');

%% tracking
uHat = zeros(parameters.simulationTime ,4);

% which trajectory estimate
trajectory = 1;

% mean of Prior, first measurement
parameters.NiterMax = 100;
[u_0,numberOfPerformedIterations] = iterativeNLS_singolarity(parameters,APhat,TYPE,R,squeeze(rho(trajectory, 1, :))');

% if we knew the velocity
x_mean = [u_0(numberOfPerformedIterations,1),u_0(numberOfPerformedIterations, 2),50, 0]';




% covariance of Prior
%P = diag(ones(1,4));
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

xlabel('[m]'), ylabel('[m]');

xlim([parameters.xmin parameters.xmax])

ylim([parameters.ymin parameters.ymax])

axis equal

grid on

% pause()





