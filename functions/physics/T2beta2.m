%====================================
% compute beta2 = v2/c2 as a function of kinetic energy
% for protons
%
% INPUT
% |T|: Kinetic energy (J)
%
% OUTPUT
% beta2=(v/c)^2
%====================================
function beta2 = T2beta2(T)

  persistent m_p c ; %Make these variables persistent between function calls to reduce time in calls to |physicsConstants|

  if ~isscalar(m_p)
    %The presistent variables have not yet been defined. Load them from physicsConstants
    physicsConstants;
  end

  %beta2 = (pc)2/E2 = (pc)2 / (T+mc2)2 % Eq 3.15 and 3.16 in [2]
  beta2 = T2pc2(T) ./ (T + m_p .* c.^2).^2; %Eq 3.19 in [2]
end
