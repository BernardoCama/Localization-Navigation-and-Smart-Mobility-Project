function plotEllipse(parameters,AP,ellipsePoints,UE,TYPE)
x_ell = ellipsePoints(1,:);
y_ell = ellipsePoints(2,:);

fig = figure(1); 
hold on
%fig.WindowState = 'maximized';
plot( AP(:,1) , AP(:,2) , '^','MarkerSize',10,'MarkerEdgeColor',[147,0,0]./255,'MarkerFaceColor',[147,0,0]./255)
fill(x_ell,y_ell,[.9 .95 1],'edgecolor',[0, 0.4470, 0.7410],'linewidth',2);alpha(.5)
%plot(UE(1),UE(2),'pk','LineWidth',1,'MarkerSize',9,'MarkerFaceColor',[.4 .4 1]);
%plot([AP(:,1),repmat(UE(1),size(AP,1),1)]',[AP(:,2),repmat(UE(2),size(AP,1),1)]','--','LineWidth',1,'color',[220,220,220]./255)
legend('AP','Ellipse','UE','location','best')
xlabel('[m]'), ylabel('[m]');
grid on
xlim([parameters.xmin parameters.xmax])
ylim([parameters.ymin parameters.ymax])
% axis equal
switch TYPE
    case 'TOA'
        title([' ',num2str(TYPE),', $N_{AP}$ = ',num2str(parameters.numberOfAP),' , $\sigma $ = ',num2str(parameters.sigmaTOA),'m'],'Interpreter','Latex')
    case 'RSS'
        title([' ',num2str(TYPE),', $N_{AP}$ = ',num2str(parameters.numberOfAP),' , $\sigma $ = ',num2str(parameters.sigmaRSS),'m'],'Interpreter','Latex')
    case 'AOA'
        title([' ',num2str(TYPE),', $N_{AP}$ = ',num2str(parameters.numberOfAP),' , $\sigma $ = ',num2str(rad2deg(parameters.sigmaAOA)),' deg '],'Interpreter','Latex')
    case 'TDOA'
        title([' ',num2str(TYPE),', $N_{AP}$ = ',num2str(parameters.numberOfAP),' , $\sigma $ = ',num2str(parameters.sigmaTDOA),'m'],'Interpreter','Latex')
end


end