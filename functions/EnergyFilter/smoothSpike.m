%% smoothSpike
% Convert the shape of a CEF spike from a step function inot a smooth function
%
%% Syntax
% |res = help_header(im1,im2)|
%
%
%% Description
% |res = help_header(im1,im2)| Description
%
%
%% Input arguments
% |im1| - _STRING_ -  Name
%
%
%% Output arguments
%
% |res| - _STRUCTURE_ -  Description
%
%
%% Contributors
% Authors : R. LAbarbe (open.reggui@gmail.com)

function p = smoothSpike(a_max , a_min , h_step)

  r = (a_max + a_min) ./2;  %Centre of the step
  h = h_step; %Height of the step
  options = optimset('Display','iter');
  %p0 = [0 , 1 , 4 , 1 , 1 , 10 , 1]
  p0 = [min(h_step) , 2 , 1 , 1 , 1];
  [p,fval,exitflag] = fminsearch(@gofSpike,p0,options, r , h)

  rsim = 0:0.1:max(r);
  f = spikeShapeR2H(p,rsim);
  figure(600)
  plot(r,h,'o')
  hold on
  plot(rsim,f,'-')
  xlabel('Step centre (mm)')
  ylabel('Step height (mm)')
  grid on
  drawnow


end


%===============================
function gof = gofSpike(p,r,y)
    f = spikeShapeR2H(p,r);
    gof = sum((y-f).^2);
end
