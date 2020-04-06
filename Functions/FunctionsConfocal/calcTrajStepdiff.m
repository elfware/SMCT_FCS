function [diffPrinc1,diffPrinc2,allDiffPrinc1,allDiffPrinc2,meanDiffPrinc1,meanDiffPrinc2]=calcTrajStepdiff(trajs,scale,deltat,inds)
% given cell of trajectoreis in pixels (first and second column positions
% in the cell array), computes the difference in the steps, returns the
% diffs both as cell array for each trajectory and as a vector with sll
% diffs for the indexes specified by inds. the diffs returned are in units
% of micrometers per second (scale is nm/ pixel and deltat time step in
% milli seconds)
numtrajs=length(trajs);
diffPrinc1=cell(1,numtrajs);
diffPrinc2=cell(1,numtrajs);
allDiffPrinc1=[];
allDiffPrinc2=[];
meanDiffPrinc1=NaN(1,numtrajs);
meanDiffPrinc2=NaN(1,numtrajs);
k=1;
for i=1:numtrajs
    
    if ~isempty(trajs{i})
    diffPrinc1{i}=diff(trajs{i}(:,1))*scale*10^-3/(deltat*10^-3);
    diffPrinc2{i}=diff(trajs{i}(:,2))*scale*10^-3/(deltat*10^-3);
    meanDiffPrinc1(i)=mean(diffPrinc1{i});
    meanDiffPrinc2(i)=mean(diffPrinc2{i});
    if k<=numel(inds)&&inds(k)==i
        allDiffPrinc1=[allDiffPrinc1;diffPrinc1{i}];
        allDiffPrinc2=[allDiffPrinc2;diffPrinc2{i}];
        k=k+1;
    end
    end
end

end

