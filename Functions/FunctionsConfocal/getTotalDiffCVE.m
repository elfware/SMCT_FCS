function [coeffDiff,coeffDiffComponents] = getTotalDiffCVE(trajs,timeStep,pixelSize,Dimension)
 %Calculculates diff coeff from trajectories according to Vestegaard et al (CVE estimator)
    %-output: coeffDiffs in um^2/s
    %trajs: cell array with where all elements have two columns with x and
    %y coordinates inp pixels
    %Timestep in milliseconds
    %pixelSize in nm
    %Dimension (2 for x and y data)
    numtrajs=length(trajs);
    coeffDiffComponents=NaN(1,Dimension);
    coeffDiff=0;
    thediffsx=[];
    thediffsy=[];
    for i=1:numtrajs
        if~isempty(trajs{i})
            thediffsx=[thediffsx;diff(trajs{i}(:,1))*pixelSize*10^-3];
            thediffsy=[thediffsy;diff(trajs{i}(:,2))*pixelSize*10^-3];
         end
        
    end
    
    
    coeffDiffComponents(1)=mean((thediffsx).^2)/(2*timeStep*10^-3)+...
        mean(thediffsx(1:end-1).*thediffsx(2:end))/(timeStep*10^-3);
    
    coeffDiffComponents(2)=mean((thediffsy).^2)/(2*timeStep*10^-3)+...
        mean(thediffsy(1:end-1).*thediffsy(2:end))/(timeStep*10^-3);
    
    coeffDiff=mean(coeffDiffComponents);
end