function [coeffDiff,coeffDiffComponents] = getDiffCVECamera(trajs,timeStep,pixelSize,Dimension)
 %Calculculates diff coeff from trajectories according to Vestegaard et al (CVE estimator)
    %-output: coeffDiffs in um^2/s
    %trajs: cell array with where all elements have two columns with x and
    %y coordinates inp pixels
    %Timestep in milliseconds
    %pixelSize in nm
    %Dimension (2 for x and y data)
numtrajs=length(trajs);
coeffDiffComponents=NaN(Dimension,numtrajs);
coeffDiff=NaN(1,numtrajs);
    for i=1:numtrajs
      if ~isempty(trajs{i})  
        for j=1:Dimension
            thediffs=diff(trajs{i}(:,j))*pixelSize*10^-3;
            coeffDiffComponents(j,i)=mean((thediffs).^2)/(2*timeStep*10^-3)+...
               mean(thediffs(1:end-1).*thediffs(2:end))/(timeStep*10^-3);
        end
        coeffDiff(i)=mean(coeffDiffComponents(:,i));
      end
    end
end