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

