%
% N36=predN36depth(pp,sp,sf,cp,age,depths,scaling_model)
%
% Predicts the cosmogenic 36-Cl concentration for a depth sample.
%
%  pp,sp,sf,cp            Chloe style physics, site, sample, and
%                         computed parameters.
%  age                    Exposure age in kyr.
%  depths                 Depths in g/cm^2.
%  scaling_model		  The scaling model being used
%
function N36=predN36depth(pp,sp,sf,cp,age,depths,erosionrate,scaling_model)
%
% Get the erosion rate in gm/(cm^2*yr) (sp.epsilon is in gm/cm^2/kyr)  
%
%erosionrate=sp.epsilon/1000;
%
% Note that from now on in this function, the erosion rate is in g/cm^2/yr.
%
%
% Adjust contemporary depth to depth at start time.
%
currentdepths=depths+erosionrate*age*1000;
%
% Pick deltat for the integration.  
%
deltat=0.1;
%
% We integrate until the point in time where the sample was collected.
%
tfinal=sp.tfinal;
%
% Figure out the depth spacing.
%
deltadepth=deltat*erosionrate*1000;
%
% Now, the main integration loop.
%
N36=0;
t=tfinal-age;
while (t+deltat < tfinal)
%
% Update the elevation/latitude scaling factor to the current
% time.  Note that our t is in kyr, with negative values in the
% past, while the TDSF stuff is done in terms of years before present. 
%
  interptime=t+deltat/2;
  sf.currentsf=getcurrentsf(sf,interptime,scaling_model,'cl');
%
% Compute the production rate.  We use the mid point of the range
% of depths corresponding to the current time interval.  
%
  pz=prodz36(currentdepths-deltadepth/2,pp,sf,cp);
%
% Update N36.
%
% There are two terms here.  The first term is the old inventory of cosmogenic
% 36-Cl, decayed for a time period of deltat.  The second term represents the
% newly generated 36-Cl, including the radioactive decay of some of it that
% was generated early in the time period.  The radioactive decay factor is
%
% f=int(P*exp(-r(deltat-t)),t=0..deltat)=P*(1-exp(-r*deltat))/r
%
% Note that for extremely small values of deltat you could get roundoff 
% errors when taking 1-exp(-r*deltat).  This hasn't been a problem for 
% our range of deltat values.
%
% The effect of using this term is that predN36's results are essentially
% indepedendent of deltat if the production rates are constant in time and
% the erosion rate is 0.   If the erosion rate is nonzero, then deltat must
% be small enough that very little erosion occurs during a time period, or
% else predN36's result will depend on deltat.  
% 
  f=(1.0-exp(-pp.lambda36Cl*deltat*1000))/pp.lambda36Cl;
  N36=N36*exp(-pp.lambda36Cl*deltat*1000)+...
      pz*f;
%
% Update t.
%
  t=t+deltat;
%
% Update depth
%
  currentdepths=currentdepths-deltat*1000*erosionrate;
end
%
% One final fractional time step.  deltatf is the length of this
% fractional time step.
%
deltatf=tfinal-t;
deltadepth=deltatf*erosionrate*1000;
%
% Update the elevation/latitude scaling factor to the current
% time.  Note that our t is in kyr, with negative values in the
% past, while the TDSF stuff is done in terms of years before present. 
%
interptime=t+deltatf/2;
sf.currentsf=getcurrentsf(sf,interptime,scaling_model,'cl');
%
% Compute the production rates.
%
pz=prodz36(currentdepths-deltadepth/2,pp,sf,cp);
%
% Update layerN36.
%
f=(1.0-exp(-pp.lambda36Cl*deltatf*1000))/pp.lambda36Cl;
N36=N36*exp(-pp.lambda36Cl*deltatf*1000)+...
    pz*f;
