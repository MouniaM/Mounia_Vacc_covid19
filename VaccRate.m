
function[VaccRate]=VaccRate(t,tVacc,DurVacc,NAge,a,AV,gamma,gamma2)
% if (t<tVacc && 1<=a<=NAge)
%       VaccRate=0;
%   if (t>=tVacc && a==AV)
    if (t>=tVacc && t<=(tVacc+DurVacc) && 1<=a<=NAge )
    VaccRate=gamma;
%     elseif (t>tVacc+365 && 1<=a<=NAge )
%     VaccRate=gamma2;
    else
        VaccRate=0;
end
return
end