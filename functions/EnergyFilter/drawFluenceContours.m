%% drawFluenceContours
% Draw the fluence contours plots
%
%% Syntax
% |drawFluenceContours(X , Y, Fluence , E_filter , figID)|
%
%
%% Description
% |drawFluenceContours(X , Y, Fluence , E_filter , figID)| Description
%
%
%% Input arguments
% |X| -_SCALAR MATRIX_- |X(x,y)| coordinate (mm) where |Fluence(x,y,E)| is computed
%
% |Y| -_SCALAR MATRIX_- |Y(x,y)| coordinate (mm) where |Fluence(x,y,E)| is computed
%
% |Fluence(x,y,E)| -_SCALAR MATRIX_-  |Fluence(x,y,E)|  fluence at position (x,y) for proton of energy |E_filter(E)|
%
% |E_filter| -_SCALAR MATRIX_- Energy (MeV) for the step(E)
%
% |figID| -_INTEGER_- ID of the first figure to draw. The other figure will increment that number. Ther is one figure per energy layer
%
%
%% Output arguments
%
% |res| - _STRUCTURE_ -  Description
%
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)

function drawFluenceContours(X , Y, Fluence , E_filter , figID , titleText)

  N_E_filter = size(Fluence,3);
  NbRows = ceil(N_E_filter./3); %3 plots per row

  figure(figID)

  for k=1:N_E_filter

      subplot(NbRows,3,k)
      if (numel(unique(Fluence(:,:,k)')) > 1)
        %Display contour only if the fluence is not uniform. Otherwise there iwll be an error message
        contour(X',Y',Fluence(:,:,k)','ShowText','on');
        colorbar
      end
      ylabel('Y (mm)')
      xlabel('X (mm)')
      title(['E = ',num2str(E_filter(k)),' MeV']);
  end
  if nargin >=6
      sgtitle(titleText)
  end

  drawnow
end
