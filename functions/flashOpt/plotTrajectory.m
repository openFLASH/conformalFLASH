%% plotTrajectory
% PLot the scanning trajectory on the BEV plot
%
%% Syntax
% |plotTrajectory(spot , beamletSequence , plotNb)|
%
%
%% Description
% |plotTrajectory(spot , beamletSequence , plotNb)| Description
%
%
%% Input arguments
% |spot| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(i,:) = [x,y]
%
% |beamletSequence| - _SCLAR MATRIX_ - The i-th spot to deliver is spot(beamletSequence(i),:)
%
% |plotNb| -_SCALAR_-  Figure number where to plot the trajectory
%
% |Tstart| -_SCALAR VECTOR_- Time stamp (ms) at the start of delivery of i-th spot
%
% |SpotIdxNoProton| -_SCLAR VECTOR_- [OPTIONAL] Index of the spots in |spot| to be displayed with a cyan *
%
% |b| -_SCALAR_- [OPTIONAL : if omited, no title is displayed] Beam number (to be used in figure title)
%
% |TimePerSpot| -_SCALAR VECTOR_- |TimePerSpot(s| Time (ms) required to deliver all the protons in s-th spot of the seqeunce. Used to represent the spot weight
%
%% Output arguments
%
% None
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function plotTrajectory(spot , beamletSequence , plotNb , Tstart , SpotIdxNoProton, b , TimePerSpot)

    if nargin < 5
      SpotIdxNoProton = [];
    end

    if nargin < 5
      b = [];
    end

    figure(plotNb)
    %plot(spot(beamletSequence(:),1),spot(beamletSequence(:),2),'ok') %Plot the spots
    minT = min(TimePerSpot(beamletSequence(:)),[],'all');
    maxT = max(TimePerSpot(beamletSequence(:)),[],'all');
    scatter(spot(beamletSequence(:),1) , spot(beamletSequence(:),2) , 50 , round(255.*(TimePerSpot(beamletSequence(:))-minT)./(maxT-minT)) , 'filled')
    hold on
    plot(spot(beamletSequence(:),1),spot(beamletSequence(:),2),'-k') %plot the trajectory line
    xlabel('X (mm)')
    ylabel('Y (mm)')
    hcb = colorbar;
    set(get(hcb,'Title'),'String','Spot charge (AU)')

    if(~isempty(b))
      title(['Beam ', num2str(b), ': Trajectory with Start Time (ms)'])
      %title(['Beam ', num2str(b), ': Trajectory with dose per spot (Gy)'])
    end

    grid on
    for idx = 1:numel(beamletSequence)
      %Indicate start time next to spot
      if(~isempty(Tstart))
        text(spot(beamletSequence(idx),1),spot(beamletSequence(idx),2), cellstr(num2str(round(Tstart(beamletSequence(idx))))))
      else
        text(spot(beamletSequence(idx),1),spot(beamletSequence(idx),2), cellstr(num2str(idx)))
      end

    end

    if(~isempty(SpotIdxNoProton))
      for idx = 1:numel(SpotIdxNoProton)
        plot(spot(SpotIdxNoProton,1),spot(SpotIdxNoProton,2),'*c') %Highlight in cyan the spots without dose
      end
end
