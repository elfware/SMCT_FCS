function [slidingVector,notstuck,angles,rellat,trajPrinclatent] = slidingClassifyPCAver(Trajectories,DataSelectionRules,classparams,p)
% Classification of trajectoreis in alltrajs, according PCA method (find
% direction of sliding as first principal compnent). Then deciding sliding
% or not based on ard cutoffs in classparams
% Change log:
% 2018-07-07: princomp is removed from Matlab 2018a and changed to pca,
% changed princomp to pca

  numtrajs=size(Trajectories.Traj,1);
  %trajIndex=1:numtrajs;
  numwindows=length(classparams.varwindow);
  minWindowVar1=zeros(numwindows,numtrajs);
  maxWindowVar1=zeros(numwindows,numtrajs);
  minWindowVar2=zeros(numwindows,numtrajs);
  maxWindowVar2=zeros(numwindows,numtrajs);
  trajL=zeros(1,numtrajs);

  angles=zeros(1,numtrajs);

  trajPrinccoeff=zeros(2,numtrajs);
  trajPrinclatent=zeros(2,numtrajs);

  %meandiffP1=zeros(1,numtrajs);
  %meandiffP2=zeros(1,numtrajs);
  %vardiffP1=zeros(1,numtrajs);
  %vardiffP2=zeros(1,numtrajs);
  %absvarP1=zeros(1,numtrajs);
  %absvarP2=zeros(1,numtrajs);
  %diffP1=cell(1,numtrajs);
  %diffP2=cell(1,numtrajs);
  thescores=cell(1,numtrajs);
  tempcoeff=cell(1,numtrajs);
  
for i=1:numtrajs
    eval(['TrajData=Trajectories.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode ';']);
    alltrajs =TrajData{i}(:,2:3);%[um]
    deltat   =diff(TrajData{i}(:,1));%[s]
    deltat(end+1)   = deltat(end);
    deltat(:,2)=deltat(:,1);

    trajL(i)=size(alltrajs,1);
    [tempcoeff{i},thescores{i},trajPrinclatent(:,i)]=pca(alltrajs./deltat);
    trajPrinccoeff(:,i)=tempcoeff{i}(:,1);
    angles(i)=angle(trajPrinccoeff(1,i)+trajPrinccoeff(2,i)*1j);
    %diffP1{i}=abs(diff(thescores{i}(:,1)));
    %diffP2{i}=abs(diff(thescores{i}(:,2)));
    %meandiffP1(i)=mean(diffP1{i});
    %meandiffP2(i)=mean(diffP2{i});
    %vardiffP1(i)=var(diffP1{i});
    %vardiffP2(i)=var(diffP2{i});
    %absvarP1(i)=var([0 cumsum(diffP1{i}')]);
    %absvarP2(i)=var([0 cumsum(diffP2{i}')]);
    
    for k=1:length(classparams.varwindow)
         window=classparams.varwindow(k);
         windowVar1=zeros(1,trajL(i)-window);
         windowVar2=zeros(1,trajL(i)-window);
         
        clear variance
        for j=1:2
	  Mean=conv(thescores{i}(:,j),ones(1,window)/window);
      Mean=Mean(window:end-window);
	  S2=conv(thescores{i}(:,j).^2,ones(1,window));
      S2=S2(window:end-window);
	  S1=conv(thescores{i}(:,j),ones(1,window));
      S1=S1(window:end-window);
	  variance(:,j)=(S2-2.*Mean.*S1+window.*Mean.^2)./(window-1);
        end
        windowVar1=variance(:,1);
        windowVar2=variance(:,2);
        
        minWindowVar1(k,i)=min(windowVar1);
        maxWindowVar1(k,i)=max(windowVar1);
        minWindowVar2(k,i)=min(windowVar2);
        maxWindowVar2(k,i)=max(windowVar2);
    end    
end

rellat=trajPrinclatent(1,:)./sum(trajPrinclatent);
notstuck=trajPrinclatent(1,:)>classparams.minNotStuck;
slidingVector=notstuck;
disp(['Not stuck: ' num2str(find(slidingVector))])

              %Select on the half circle -pi/2 to pi/2                                          %Change from the interval [-pi pi] to [0 2*pi] and then select on the interval pi/2 to 3pi/2
AngleSelected=((angles>classparams.AngleInterval(1)) & (angles<classparams.AngleInterval(2))) | ((mod(angles,2*pi)>mod(classparams.AngleInterval(1)+pi,2*pi)) & (mod(angles,2*pi)<mod(classparams.AngleInterval(2)+pi,2*pi)));
slidingVector=slidingVector & AngleSelected;
disp(['Angle: ' num2str(find(slidingVector))])

rellatvector=rellat>classparams.PrincRelVar(1) & rellat<classparams.PrincRelVar(2);
slidingVector=slidingVector & rellatvector;
disp(['PrincRelVar: ' num2str(find(slidingVector))])

trajPrinclatentVector=trajPrinclatent(1,:)>classparams.PrincAbsVar(1) & trajPrinclatent(1,:)<classparams.PrincAbsVar(2);
slidingVector=slidingVector & trajPrinclatentVector;
disp(['PrincAbsVar: ' num2str(find(slidingVector))])

%slidingInd=slidingInd(minWindowVar1(1,slidingInd)>classparams.minMinWindowVar);
%disp(['MinMinWindowVar: ' num2str(slidingInd)])

if p
  figure()
  subplot(1,2,1)
  plot(angles(notstuck & not(slidingVector)),rellat(notstuck & not(slidingVector) ),'bx');
  hold on;plot(angles(slidingVector),rellat(slidingVector),'rx')
  plot([-4 4],classparams.PrincRelVar(1)*[1 1],'-r')
  plot([classparams.AngleInterval(1) classparams.AngleInterval(1)],[0 max(rellat(notstuck))],'-b')
  plot([classparams.AngleInterval(2) classparams.AngleInterval(2)],[0 max(rellat(notstuck))],'-b')
  plot([classparams.AngleInterval(1)+pi classparams.AngleInterval(1)+pi],[0 max(rellat(notstuck))],'-b')
  plot([classparams.AngleInterval(2)-pi classparams.AngleInterval(2)-pi],[0 max(rellat(notstuck))],'-b')
  
  Values=find(notstuck);
  for I=1:size(Values,2) 
    text(angles(Values(I))+0.001,rellat(Values(I)),num2str(Values(I)))
  end
  axis([-4 4 0 max(rellat(notstuck))])
  xlabel('Angle of princ componenet with largest eig');
  ylabel('Relative variance in Princ direction');
  title(['Tot num of traj is ' num2str(size(notstuck,2))])
  hold off
  
  subplot(1,2,2)
  plot(angles(notstuck & not(slidingVector)),trajPrinclatent(1,notstuck & not(slidingVector)),'bx')
  hold on;plot(angles(slidingVector),trajPrinclatent(1,slidingVector),'rx')
  plot([-4 4],classparams.PrincAbsVar(1)*[1 1],'-r')
  plot([classparams.AngleInterval(1) classparams.AngleInterval(1)],[0 max(trajPrinclatent(1,notstuck))],'-b')
  plot([classparams.AngleInterval(2) classparams.AngleInterval(2)],[0 max(trajPrinclatent(1,notstuck))],'-b')
  plot([classparams.AngleInterval(1)+pi classparams.AngleInterval(1)+pi],[0 max(trajPrinclatent(1,notstuck))],'-b')
  plot([classparams.AngleInterval(2)-pi classparams.AngleInterval(2)-pi],[0 max(trajPrinclatent(1,notstuck))],'-b')

  Values=find(notstuck);
  for I=1:size(Values,2) 
    text(angles(Values(I))+0.001,trajPrinclatent(1,Values(I)),num2str(Values(I)))
  end
  axis([-4 4 0 max(trajPrinclatent(1,notstuck))])
  xlabel('Angle of princ componenet with largest eig');
  ylabel('variance in Princ direction');
  hold off
  
  k=0;
  if k
    index=9; %5 6 19 20 33 56
    figure()
    scatter(Trajectories.TrajAvriage{index}(:,2) , Trajectories.TrajAvriage{index}(:,3))
    figure()
    scatter(thescores{index}(:,1) , thescores{index}(:,2))
    %ShowTraj(Scanobject,Trajectories,DataSelectionRules, index)
    tempcoeff{index}
    trajPrinclatent(:,index)
    angles(index)
  end
end
