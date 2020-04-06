function [coeffDiff,numberOfCells] = getMSD(outputFolder,timeStep,pixelSize,numberOfFittingPoints,Dimension)
% This function computes apparent diffusion coefficients for trajectories
% of all sizes
%
% Input:
%   - outputFolder : folder that contains the results
%   - timeStep : time between frames [in ms]
%   - pixelSize : size of the pixels [in nm] (typical value is 80)
%   - numberOfFittingPoints : number of points to fit to get the slope
%   - Dimension : 2 for 2D trajectories, 1 for 1D trajectories, etc.
%
% Output:
%   - list of diffusion coefficients
%   - plot histogram for diffusion coefficients
%
% Author : Vladimir Curic, 2015-01-05

positionList = dir(strcat(outputFolder,filesep,'Pos*'));
coeffDiff = [];
numberOfCells = 0;
for positionIdx = 1:length(positionList)
    disp(positionList(positionIdx).name)
    
    cellMatList = dir(strcat(outputFolder,filesep,positionList(positionIdx).name,filesep,'*.mat'));
    cellMatFile = load(strcat(outputFolder,filesep,positionList(positionIdx).name,filesep,cellMatList(1).name));
    fieldName = fieldnames(cellMatFile);
    cellStruct = cellMatFile.(fieldName{1});
    
    numberOfCells = numberOfCells + size(cellStruct,2);
    
    for jj = 1:length(cellStruct)
        
        numberTempTrajectories = size(cellStruct(jj).particles,2);
        meanDispl = [];
        for kk = 1:numberTempTrajectories
            tempVar = cellStruct(jj).particles;
            tempTraj = tempVar(kk).trajImageCoord.*0.08;
            nData = size(tempTraj,1); %# number of data points
            numberOfdeltaT = nData-1;
            msd = [];
            for dt = 1:numberOfdeltaT
                deltaCoords = tempTraj(1+dt:end,1:2) - tempTraj(1:end-dt,1:2);
                squaredDisplacement = sum(deltaCoords.^2,2); %# dx^2+dy^2+dz^2
                
                msdMean = mean(squaredDisplacement); %# average
                msdStd = std(squaredDisplacement); %# std
                msdN = nData-length(squaredDisplacement); %# n
                msd = [msd; msdMean, msdStd, msdN, kk];
            end
            
            tempFit = polyfit(msd(1:numberOfFittingPoints,3),msd(1:numberOfFittingPoints,1),1);
            
            coeffDiff = [coeffDiff tempFit(1,1)/(2*Dimension*timeStep*10^(-3))];        
            
        end
    end
end

figure, hist(coeffDiff,50)
%title('Diffusion coefficients')
xlabel('Apparent diffusion coefficient [$\mu m^2 s^{-1}$]','Interpreter','latex')
ylabel('Number of of trajectories','Interpreter','latex')

end



