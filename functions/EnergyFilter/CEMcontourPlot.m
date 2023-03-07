%% CEMcontourPlot
% Draw a contour plot of the elevation map of the CEM.
% Optionally overlay the aperture contour on the elevation map. The aperture contour is projected onto the CEM plane.
%
%% Syntax
% |CEMcontourPlot(fig , Xvec, Yvec, CEFelevation , BlockData)|
%
%
%% Description
% |CEMcontourPlot(fig , Xvec, Yvec, CEFelevation )| Display elevation map of the CEM
%
% |CEMcontourPlot(fig , Xvec, Yvec, CEFelevation , BlockData , VDSA , DisoPlane)| Display elevation map of the CEM and the aperture contour
%
%% Input arguments
% |fig| - _SCALAR_ - Matlab figure number
%
% |Xvec| -_SCALAR VECTOR_- x IEC gantry coordinate of the point at |CEFelevation(x,:)|
%
% |Yvec| -_SCALAR VECTOR_-  y IEC gantry coordinate of the point at |CEFelevation(:,j)|
%
% |CEFelevation| -_SCALAR MATRIX_- |CEFelevation(x,y)| Thickness (mm) of the CEF pixel at position (x,y) in the IEC beam Limiting device CS
%
% |BlockData| -_SCALAR MATRIX_- [OPTIONAL. If absent, aperture is not displayed] |BlockData{cntIdx}(i,:)=[x,y]|  Coordinates (mm) of the i-th point defining the contour of the aperture block projected onto the machine isocentric plane in the IEC BEAM LIMITING DEVICE coordinate system for the b-th beam.
%
% |VDSA| -_SCALAR VECTOR_- [OPTIONAL. Only required if |BlockData| is present] [sadX , sadY]  mm source to axis distance for the 2 scanning magnets
%
% |IsocenterToRangeModulatorDistance| -_SCALAR_- [OPTIONAL. Only required if |BlockData| is present] Distance (mm) from isocentre to base of CEM
%
%
%% Output arguments
%
% None
%
%
%% Contributors
% Authors : R. Labarbe, Lucian Hotoiu (open.reggui@gmail.com)

function CEMcontourPlot(fig , Xvec, Yvec, CEFelevation , BlockData , VDSA , IsocenterToRangeModulatorDistance)

  if nargin < 5
    BlockData = [];
  end


  figure(fig)
  fprintf('Preparing contour plot \n')
  imagesc(Xvec , Yvec , CEFelevation'); %imagesc plots the first index of |CEFelevation| vertically, from top to bottom of image
  set (gca,'Ydir','normal') %Set vertical axis in the upward direction
  hold on
  colormap Turbo
  grid on

  if ~isempty(BlockData)

    mag = project2Plane(VDSA , IsocenterToRangeModulatorDistance);

    for cntIdx = 1:numel(BlockData)
      plot(BlockData{cntIdx}(:,1).* mag(1) , BlockData{cntIdx}(:,2) .* mag(2) ,'-r' ) %imagesc plots the first index of |CEFelevation| vertically, from top to bottom of image
      hold on
    end
  end

  title(['Conformal Energy Filter @ CEF plane - iso view'])
  xlabel('X_g (mm)') %imagesc plots the first index of |CEFelevation| vertically, from top to bottom of image
  ylabel('Y_g (mm)')
  zlabel('Height (mm)')
  hcb = colorbar;
  set(get(hcb,'Title'),'String','Height (mm)')
  grid on
  drawnow

%   figure(fig+1)
%   surf(Xvec, Yvec, CEFelevation');
%   set (gca,'Ydir','normal') %Set vertical axis in the upward direction
%   hold on
%   title(['Conformal Energy Filter 3D - iso view'])
%   xlabel('X_g (mm)') %imagesc plots the first index of |CEFelevation| vertically, from top to bottom of image
%   ylabel('Y_g (mm)')
%   zlabel('Height (mm)')
%   hcb = colorbar;
%   set(get(hcb,'Title'),'String','Height (mm)')
%   grid on
%   drawnow

  figure(fig+2)
  tmp = flip(CEFelevation',1);
  tmp = flip(tmp,1);
  h = bar3(Yvec , tmp , 1 , 'y');
  %Set the correct scale for X axis
  Xdat=get(h,'XData');
  axis tight
  xstp = min(diff(Xvec));
  minX = min(Xvec);
  for ii=1:length(Xdat)
      Xdat2 = (Xdat{ii}-1) .* xstp + minX .* ones(size(Xdat{ii})) ;
      set(h(ii),'XData',Xdat2);
  end

  set (gca,'Ydir','normal') %Set vertical axis in the upward direction
  hold on
  title(['Conformal Energy Filter 3D - iso view'])
  xlabel('X_g (mm)') %imagesc plots the first index of |CEFelevation| vertically, from top to bottom of image
  ylabel('Y_g (mm)')
  zlabel('Height (mm)')
  grid on
  drawnow

end
