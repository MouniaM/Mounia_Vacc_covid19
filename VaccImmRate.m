function[VaccImmRate]=VaccImmRate(t,tVacc,NAge,a,AV,omega)
% if (t<tVacc && 1<=a<=NAge)
%       VaccImmRate=0;
%  if (t>=tVacc && a==AV)
   if (t>=tVacc && 1<=a<=NAge )
     VaccImmRate=omega;
else
    VaccImmRate=0;
end
    return
end