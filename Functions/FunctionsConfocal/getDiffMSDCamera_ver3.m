function [msd,coeffDiffComponents,coeffDiff] = getDiffMSDCamera_ver2(trajs,timeStep,pixelSize,Dimension,startPoint,endPoint)
 %Calculculates diff coeff from trajectories according to MSD curve fitting
    %-output: coeffDiffs in um^2/s
    %trajs: cell array with where all elements have two columns with x and
    %y coordinates inp pixels
    %Timestep in milliseconds
    %pixelSize in nm
    %Dimension (2 for x and y data)
numtrajs=length(trajs);
coeffDiffComponents=zeros(Dimension,numtrajs);
coeffDiff=zeros(1,numtrajs);
msd=cell(1,Dimension);
for j=1:Dimension
    msd{j}=zeros(numtrajs,endPoint);
end
    for i=1:numtrajs
        for j=1:Dimension
%             thediffs=diff(trajs{i}(:,j))*pixelSize*10^-3;
%             coeffDiffComponents(j,i)=mean((thediffs).^2)/(2*timeStep*10^-3)+...
%                mean(thediffs(1:end-1).*thediffs(2:end))/(timeStep*10^-3);
%            
           nData = size(trajs{i}(:,j),1); %# number of data points
 %          numberOfdeltaT = nData-1;
           tempTraj=trajs{i}*pixelSize*10^-3;
           for dt = 1:endPoint
               
                   deltaCoords = tempTraj(1+dt:end,j) - tempTraj(1:end-dt,j);
               

              
               squaredDisplacement = sum(deltaCoords.^2,2); %# dx^2+dy^2+dz^2
              
                msdMean = mean(squaredDisplacement); %# average
%                msdStd = std(squaredDisplacement); %# std
%                msdN = nData-length(squaredDisplacement); %# n
               msd{j}(i,dt) = msdMean;
           end
            t=startPoint*timeStep:timeStep:endPoint*timeStep;
            tempfit=polyfit(t/1000,msd{j}(i,startPoint:endPoint),1);
            coeffDiffComponents(j,i)=tempfit(1);
        end
        
         
         
           for dt = 1:endPoint
               
                   deltaCoords = tempTraj(1+dt:end,1:2) - tempTraj(1:end-dt,1:2);
               

              
               squaredDisplacement = sum(deltaCoords.^2,2); %# dx^2+dy^2+dz^2
              
                msdMean = mean(squaredDisplacement); %# average
%                msdStd = std(squaredDisplacement); %# std
%                msdN = nData-length(squaredDisplacement); %# n
               msd{3}(i,dt) = msdMean;
           end
            t=startPoint*timeStep:timeStep:endPoint*timeStep;
            tempfit=polyfit(t/1000,msd{j}(i,startPoint:endPoint),1);
            coeffDiff(i)=tempfit(1);
        
%         coeffDiff(i)=mean(coeffDiffComponents(:,i));
    end
end