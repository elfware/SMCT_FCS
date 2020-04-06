function Trajectories=getDiffCVETotal(Trajectories,DataSelectionRules,p)
 %Calculculates diff coeff from trajectories according to Vestegaard et al (CVE estimator)
    %-output: coeffDiffs in um^2/s
    %trajs: cell array with where all elements have two columns with x and
    %y coordinates in micrometers
    %Timestep in seconds
    
%numtrajs=length(trajs);

%fieldnames(Trajectories.TrajPosEstimation){fieldIndex}
numtrajs=size(Trajectories.Traj,1);
    for i=1:numtrajs
        eval(['Data=Trajectories.TrajPosEstimation.' DataSelectionRules.Trajectory.Mode ';']) 
        thediffs=diff(Data{i}(:,2:3),1,1);%[um]
        timeStep=diff(Data{i}(:,1),1,1);%[s]
            
        coeffDiffComponents(i,:)= mean(thediffs.^2./(2*repmat(timeStep,1,2))) + mean(thediffs(1:end-1,:).*thediffs(2:end,:)./repmat(timeStep(1:end-1),1,2));   
    end
%coeffDiff=mean(coeffDiffComponents,2); Hmm this can not be right? 
coeffDiff=sqrt(coeffDiffComponents(:,1).^2 + coeffDiffComponents(:,2).^2);

Trajectories.Traj(:,end+[1:2])=coeffDiffComponents;
Trajectories.Traj(:,end+1)=coeffDiff;

S=size(Trajectories.Traj,2);

Trajectories.INFO.Traj=[Trajectories.INFO.Traj;...
			[num2str(S-2) '.Diffusion coefficent x [um^2/s]'];...
			[num2str(S-1) '.Diffusion coefficent y [um^2/s]'];...
			[num2str(S)   '.Diffusion coefficent sqrt(x^2+y^2) [um^2/s]']];


if p>0
  coeffDiffx=coeffDiffComponents(:,1);
  coeffDiffy=coeffDiffComponents(:,2);

  if p==1
  figure()
  for Index=[1:size(Trajectories.INFO.Traj,1)]     
      if not(isempty(strfind(Trajectories.INFO.Traj{Index},'Mean polarization')))
          dataplace=str2double(Trajectories.INFO.Traj{Index}(1));
      end
  end
  
  scatter(coeffDiffx(:),coeffDiffy(:),5*ones(size(coeffDiffy(:))),Trajectories.Traj(:,dataplace),'filled')
  %plot(coeffDiffx(:),coeffDiffy(:),'.');
  title('Diffution constant map')
  xlabel('Dx [um^2/s]');
  ylabel('Dy [um^2/s');
  grid on
  axis([min([coeffDiffx(:);coeffDiffy(:)]) max([coeffDiffx(:);coeffDiffy(:)]) min([coeffDiffx(:);coeffDiffy(:)]) max([coeffDiffx(:);coeffDiffy(:)])])
  colorbar
  hold on 
  
  for I=1:size(coeffDiffx,1) 
    text(coeffDiffx(I)+0.001,coeffDiffy(I),num2str(I))
  end
  hold off
  end
  
  if p==2
  figure()
  for Index=[1:size(Trajectories.INFO.Traj,1)]
      if not(isempty(strfind(Trajectories.INFO.Traj{Index},'Mean polarization')))
          dataplace=str2double(Trajectories.INFO.Traj{Index}(1));
      end
  end
  scatter(abs(coeffDiffx(:)),abs(coeffDiffy(:)),5*ones(size(coeffDiffy(:))),Trajectories.Traj(:,dataplace),'filled')
  %plot(coeffDiffx(:),coeffDiffy(:),'.');
  title('abs Diffution constant map')
  xlabel('Dx [um^2/s]');
  ylabel('Dy [um^2/s');
  grid on
  axis([0 max(abs([coeffDiffx(:);coeffDiffy(:)])) 0 max(abs([coeffDiffx(:);coeffDiffy(:)]))]);
  colorbar
  hold on 
  
  for I=1:size(coeffDiffx,1) 
    text(abs(coeffDiffx(I))+0.001,abs(coeffDiffy(I)),num2str(I));
  end
  hold off
  end
end
