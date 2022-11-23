%% EvsDepth
% Compute the remaining energy (MeV) of the proton that has travelled to depth X (cm) in a material
%
%% Syntax
% |Er = EvsDepth(x,alpha,p, E)|
%
%
%% Description
% |Er = EvsDepth(x,alpha,p, E)| Description
%
%
%% Input arguments
%
%  |x| - _SCALAR VECTOR_ -  depth in the material (cm)
%  |alpha| - _SCALAR_ - (cm MeV^(-p)) Factor of the  range vs energy curve (Bragg Kelman equation).
%  |p| - _SCALAR_ - Exponent of the Bragg Kelman equation. No unit.
%  |E| - _SCALAR_ - Kinetic energy (MeV) of the incoming particle
%
%% Output arguments
%
%  |Er|- _SCALAR VECTOR_ - |Er(i)| Kinetic energy (MeV) of the particle at depth |x(i)| inside the material
%
%% Contributors
% Authors : R. Labarbe (open.reggui@gmail.com)
%
% REFERENCE
% [1] Bortfeld, T. (1997). An analytical approximation of the Bragg curve for therapeutic proton beams, 2024–2033. Retrieved from https://physics.fjfi.cvut.cz/~contrgui/IJZ/Lectures_2017/Bortfeld.pdf


function Er = EvsDepth(x,alpha,p, E)
  R0 = energy2range(E, alpha,p); %Maximum range of the incoming particle (cm)
  Er = range2energy(R0-x, alpha,p); %Kinetic energy (MeV) remaining at depth x

  %Negative range means that the particle stopped inside the material
  wneg = find((R0-x) < 0);
  Er(wneg) = 0;
end
