function [ h_uhat ] = createVectorOfObservations(parameters,uHat,AP,TYPE)

switch TYPE
      case 'TOA'
          
            for i = 1:parameters.numberOfAP
  
                  h_uhat(i) = sqrt(sum([uHat-AP(i,:)].^2,2));
                  
            end %i

      case 'AOA'
          
            for i = 1:parameters.numberOfAP
                
                 h_uhat(i) = atan2((uHat(2)-AP(i,2)),(uHat(1)-AP(i,1))); 
                  
            end %i   
end 

end

