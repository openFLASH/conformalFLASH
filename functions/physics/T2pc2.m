%============================
% Compute (pc) from kinetic energy using relativistic formula
% for proton
% Input
% |T| Kinetic energy (J)
%
% OUPUT
% |pc| (J)
%============================
function pc2 = T2pc2(T)
  % (pc)2 + (mc2)2 = E2 = (T + mc2)2 % From eq 3.16 and 3.17 in [2]
  % (pc)2 + (mc2)2 = T2 + (mc2)2 + 2 T mc2
  % (pc)2  = T2 + 2 T mc2
  persistent m_p c ; %Make these variables persistent between function calls to reduce time in calls to |physicsConstants|

  if ~isscalar(m_p)
    %The presistent variables have not yet been defined. Load them from physicsConstants
    physicsConstants;
  end

  pc2  = T.^2 + 2.* T .* m_p .* c.^2;
end
