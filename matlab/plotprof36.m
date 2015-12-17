function plotprof36(comp,uncerts,depths,er,age,inher,scaling_model)
%
%   plotprof36(comp,uncerts,depths,er,age,inher)
%
%   Plots our model at given parameters vs. data

numbersamples=size(comp,1);
deepest=max(depths);
%deal with both positive and negative erosion rates

if er>=0;
    mydepths=(0.0:10:deepest)';
else
    mydepths=(-er*age:10:deepest)';
end
    
%compute average composition
averagesample=comp(1,:);
for i=13:37;
averagesample(i)=mean(comp(:,i));
end

maxdepth36 = 1000+deepest+er*age;

%compute predicted production at depth for average composition
[pp,sp36,sf36,cp36]=getpars36(averagesample,maxdepth36);
sp36.epsilon=er;

size(mydepths)
size(er)
size(age)

[predictedN36]=predN36depth(pp,sp36,sf36,cp36,age,mydepths,er/1000,scaling_model);
%[predN36measure]=predN36depth(pp,sp36,sf36,cp36,age,depths);
Prodtotal=prodz36(0,pp,sf36,cp36);
predictedN36=predictedN36+inher*Prodtotal;
%predN36measure=predN36measure+inher*Prodtotal;

%compute predicted production at each actual sample depth (prodz)
prodzpred=prodz36(depths,pp,sf36,cp36);

%compute prodz for each actual composition & actual depth
for j=1:numbersamples;
    sampledepth=depths(j);
    [pp,sp36,sf36,cp36]=getpars36(comp(j,:),(sampledepth+200));
    prodzactual(j)=prodz36(sampledepth,pp,sf36,cp36);
    %calculate normalization factor
    normalization(j)=prodzactual(j)/prodzpred(j);
end

%normalize each measured N36 by the normalization factor
%measured=zeros(1,numbersamples);
measured=comp(:,1)';
N36plot=measured./normalization;
uncert=uncerts(:,1);
%plot each normalized N36

%plot profile

%conc=comp(:,1);
%maxdepth36 = 1000+deepest+er*age;
%[pp,sp36,sf36,cp36]=getpars36(comp(1,:),maxdepth36);
%sp36.epsilon=er;

% Calculate theoretical profile
%[predictedN36]=predN36depth(pp,sp36,sf36,cp36,age, ...
%					   mydepths,scaling_model);
%Prodtotal=prodz36(0,pp,sf36,cp36);
%predictedN36=predictedN36+inher*Prodtotal;


%  Plot model with data
set(axes,'XDir','default', 'YDir', 'reverse');
% set(gca,'XAxisLocation','top')
hold on
plot(N36plot,depths,'ko');
herrorbar(N36plot,depths,uncert,'ko');
hConc = plot(N36plot,depths,'bo');
hPred = plot(predictedN36,mydepths,'k');
% set(hConc,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[.7 .7 .7])
% hTitle=title('Fitted Depth Profile');
hLegend=legend([hPred,hConc],'Predicted ^{36}Cl','Measured ^{36}Cl','Location','SouthEast');
hXLabel=xlabel('Concentration of ^{36}Cl');
hYLabel=ylabel('Depth (g/cm^2)');
plotstyle(hXLabel,hYLabel,hLegend);
hold off
