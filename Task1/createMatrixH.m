%{
function [H] = createMatrixH(parameters,uHat,AP,TYPE)
switch TYPE
      case 'TOA'
          
            for i = 1:parameters.numberOfAP
                
                  for j = 1:length(uHat)   

                        H(i,j) = -(AP(i,j)-uHat(j))/sqrt(sum([uHat-AP(i,:)].^2,2));
                  
                  end %j
                  
             end %i

      case 'AOA'
          
            for i = 1:parameters.numberOfAP
                
                  for j = 1:length(uHat)   

                        if j ==1 % x coordinate
                            H(i,j) = (AP(i,2)-uHat(2))/(sum([uHat-AP(i,:)].^2,2));

                        else     % y coordinate 
                            H(i,j) = -(AP(i,1)-uHat(1))/(sum([uHat-AP(i,:)].^2,2));

                        end

                  end %j
                  
            end %i   
end
end
%}

function [H] = createMatrixH(parameters,UE,AP,TYPE)

distanceUEAP = sqrt(sum([UE-AP].^2,2));



directionCosineX = (UE(1)-AP(:,1))./ distanceUEAP;
directionCosineY = (UE(2)-AP(:,2))./ distanceUEAP;

H= zeros(parameters.numberOfAP,2);

for a= 1: parameters.numberOfAP

      switch TYPE
            case 'TOA'
                  H(a,:) = [directionCosineX(a), directionCosineY(a)];

            case 'AOA'
                  H(a,:) = [ - directionCosineY(a)/ distanceUEAP(a), 
                          directionCosineX(a)/ distanceUEAP(a)];
      end

end

