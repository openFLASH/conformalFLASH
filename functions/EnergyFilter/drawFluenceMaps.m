%% drawFluenceMaps
% Draw the fluence map
%
%% Syntax
% |drawFluenceMaps(X , Y, Fluence , E_filter , figID)|
%
%
%% Description
% |drawFluenceMaps(X , Y, Fluence , E_filter , figID)| Description
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

function drawFluenceMaps(X , Y, Fluence , E_filter , figID , titleText)

  N_E_filter = size(Fluence,3);
  NbRows = ceil(N_E_filter./3); %3 plots per row

  figure(figID)

  for k=1:N_E_filter

      subplot(NbRows,3,k)
      imagesc(X,Y,Fluence(:,:,k)');
      %imagesc takes first index as Y and second index as X, so transpose
      ylabel('Y (mm)')
      xlabel('X (mm)')
      title(['E = ',num2str(E_filter(k)),' MeV']);
      grid on
      colorbar
  end
  if nargin >=6
      sgtitle(titleText)
  end

  drawnow

end
