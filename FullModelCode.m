function out=YWdefsolve7()

%% Experiment parameters

E=2.1*1e6;% PDMS moudlus---unit: Pa
hinitial=;% Initial central separation---unit:m
eta= 0.5;% Viscosity---unit: pa.s
R=0.0071;% Radius of curvature---unit:m
k=1021; % spring constant---unit:N/m

v=10E-6;% %% Velocity Sign: 1) postive: retraction; 2) negative: approach
poi=0.5; % Poisson's ratio

ShiftLength = 0*1e-9;
Hamaker = 10e-20;
electrolyte_conc = 0.001*6.022e23*1000;  %1mM in 1/m^3
surface_potential = -0.05;              %V
boltzman_roomtemp = 4.14e-21;           %J
electron_charge = -1.602e-19;           %C
permittivity = 8.854e-12;               %C^2/J/m
dielectric = 80;                        %unitless
debye_length =(permittivity*dielectric*boltzman_roomtemp/(2*electrolyte_conc*electron_charge^2))^0.5;
w_doublelayer = 64*electrolyte_conc*boltzman_roomtemp*tanh(electron_charge*surface_potential/(4*boltzman_roomtemp));


%% Computation domain

Rs=[0:R/5000:R/10]; % mesh size in radius
hguess1=zeros(1,length(Rs));% mesh size for hguess
p1guess=zeros(1,length(Rs));%mesh size for pressure p
hcalc=zeros(1,length(Rs));% mesh size for hcalc
wguess=zeros(1,length(Rs));% mesh size for deformation w
disjoining_p_VDW = zeros(1,length(Rs));% mesh size for Van der waals force
disjoining_p_DL = zeros(1,length(Rs));% mesh size for double layer force
allruns = zeros(400,1);%matrix to store run number
allforce=zeros(400,1);%matrix to store force 
allhvst=zeros(400,length(Rs)); %matrix to store the central separation h
allwvst=zeros(400,length(Rs)); %matrix to store the deformation w
allpvst=zeros(400,length(Rs));%matrix to store the pressure
allxvst=zeros(400,length(Rs));%matrix to store x=h-w
allVDW=zeros(400,length(Rs));%matrix to store Van der Waals force
allDL=zeros(400,length(Rs));%matrix to store double layer force
allhcalc = zeros(3000,length(Rs));%matrix to store hcalc
allhguess1 = zeros(3000,length(Rs));%matrix to store hguess
%% Non-dimensionlization
YvsR(1,:) = Rs.^2/(2*R) + hinitial;% h in radial profile
Rsdless = Rs/((R*hinitial)^0.5);% non-dimesionlized radius
YvsRdless = YvsR/hinitial;%non-dimesnialized initial separation
ShiftLengthdless = ShiftLength/hinitial;%non-dimesnialized shiftlength
thicknessmax = R;
thicknessdless = thicknessmax/((R*hinitial)^0.5);%non-dimesnialized film thickness
time=0;
tincre = ;% time increment
tincredless = tincre*-v/hinitial;%non-dimesnialized time step
time=time+tincredless; 
timestop=100;
springparam = k*(hinitial^2)/(-v*eta*(R^2));
elasticparam = eta*(R^1.5)*-v/((E/(1-poi^2))*(hinitial^2.5));
%% Initialization 
h0=YvsRdless(1,:);
allhvst(1,:) = YvsRdless(1,:);
allxvst(1,:) =YvsRdless(1,:);
hguess1=YvsRdless(1,:);
hcalc=YvsRdless(1,:);
hpreviousSol=YvsRdless(1,:);




ratio = 0.965;
alldata = [time timestop hinitial ratio eta Hamaker electrolyte_conc surface_potential tincre thicknessmax R v k poi E];
%% define Hankel transform variable
inc = 100;
maxlength = 100000;
range = 1:inc:maxlength;
rangedless = range*((R*hinitial)^0.5); %dimensionless hankel transform varialbe
Z = zeros(1,maxlength/inc);


VDW = false;
DL  = false;



%% main code
%% Iterate parameter
indexer=2;

while time>-650 %The outer loop, when the loop execute once one certain time separation was determined
    if (hcalc(1))*hinitial < 1e-9
        timestop = min(timestop, time);
        disp('stop');
    end

  
  
    hpreviousSol=hcalc;
    criteriacount = 0;
    if indexer > 2 % This if statement is for the updating initialize dh.
        hguess1=2*hcalc - allhvst(indexer-2,:);
    else
        hguess1=2*hcalc - allhvst(indexer-1,:);
    end
    
    
    criteria=false;
    
    runs=1;
    while criteria==false  % The inner loop
        %% calculate hydrodynamic pressure from hguess and lubrication equation
        dhdt1=(hguess1-hpreviousSol)./tincredless;
        inside = dhdt1.*Rsdless;
        fittinginside = spline(Rsdless,inside);      
        integ = fnint(fittinginside);
        calcint = fnval(integ,Rsdless);

        for a = 1:length(p1guess)
            dp1finder(a) = -12*calcint(a)/(Rsdless(a)*(hguess1(a))^3);
            dp1finder(1) = 0;
        end
           p1guess(length(Rs)) = (3*(allxvst(indexer-1,1)-(hguess1(1)-wguess(1)*elasticparam))/tincredless)/(hcalc(1)-wguess(1)*elasticparam+ max(Rsdless)^2/2)^2;

        for a=length(p1guess)-1:-1:1
            p1guess(a)=p1guess(a+1)+((dp1finder(a)+dp1finder(a+1))/2)*(Rsdless(a+1)-Rsdless(a));%edge is zero, matrix of dp converted to pressure profile
                  
        end
        p1guess(1) = p1guess(2);

        if VDW == true
            disjoining_p_VDW = Hamaker./(6*pi)./((((hguess1+hpreviousSol)*hinitial)*0.5).^3);
            disjoining_p_VDW = disjoining_p_VDW/(eta*R*-v/(hinitial^2));
        end
        
        if DL == true
            disjoining_p_DL = w_doublelayer*exp(-((hguess1+hpreviousSol)*hinitial)/(debye_length*2));
            disjoining_p_DL = disjoining_p_DL/(eta*R*-v/(hinitial^2));
        end
      
        
        %% calculate hcalc according to force balance and verify the difference is within the tolerance

        for i = 1:length(rangedless)
            Z(1,i) = trapz(Rsdless,Rsdless.*p1guess.*besselj(0,rangedless(1,i)*Rsdless));
        end
        
        gamma = 3-4*poi;
        X = (gamma*(1-exp(-4*rangedless*thicknessdless))-4*rangedless.*thicknessdless.*exp(-2*rangedless*thicknessdless))./(gamma*(1+exp(-4*rangedless*thicknessdless))+(gamma^2+1+4*(rangedless*thicknessdless).^2).*exp(-2*rangedless*thicknessdless));
        
        
        for j = 1:length(Rsdless)
            wguess(1,j) = trapz(rangedless,2*(1-poi.^2).*X.*Z.*besselj(0, rangedless.*Rsdless(1,j)));
          
        end
        force=trapz(Rsdless,2*pi*Rsdless.*(p1guess+disjoining_p_DL-disjoining_p_VDW))+(allxvst(indexer-1,1)-(hguess1(1)-wguess(1)*elasticparam))/(hcalc(1)-wguess(1)*elasticparam+max(Rsdless)^2/2);
        hcalc=force/springparam-min(time,timestop)+h0+elasticparam*wguess;
        xxx = hcalc-wguess*elasticparam;

        
        
        if max(abs(hcalc-hguess1))*hinitial< 1E-10
            criteria=true;
            allruns(indexer,1) = runs;
            
            
        else
            
            allhcalc(runs,:) = hcalc;
            allhguess1(runs,:) = hguess1;
      %% A method to update new hguess to decrease the difference between hcalc and hguess      
            hguess1=hguess1*ratio+hcalc*(1-ratio);
          
            if min(hcalc)<0 || min(hguess1)<0
                 ratio = ratio + (1-ratio)*0.0008;
                 hguess1 = allhguess1(runs,:)*ratio+allhcalc(runs,:)*(1-ratio);
                 criteriacount = 0;
             end
 
             if runs > 1
                  
                    if (min(hguess1-allhguess1(runs-1)) > 0) && (min(allhguess1(runs-1,:)-allhguess1(runs,:)) > 0)
                        ratio = ratio + (1-ratio)*0.0008;
                        hguess1=allhguess1(runs,:)*ratio+allhcalc(runs,:)*(1-ratio);
                        criteriacount = 0;
                    elseif (max(hguess1-allhguess1(runs-1,:)) < 0) && (max(allhguess1(runs-1,:)-allhguess1(runs,:)) < 0)
                        ratio = ratio + (1-ratio)*0.0008;
                        hguess1=allhguess1(runs,:)*ratio+allhcalc(runs,:)*(1-ratio);
                        criteriacount = 0;
                    else
                        criteriacount = criteriacount +1;
                    end
                    

                    if (criteriacount >350) && (max(abs(hcalc-hguess1)*hinitial) <1E-8)
                        ratio = ratio-ratio/20000;
                        criteriacount = 0;
                    end

             end
%% plot of the results
            subplot(1,4,3)
            plot(Rs,(allhcalc(runs,:)-allhguess1(runs,:))*hinitial)
            title(num2str(ratio));
             axis([0 0.002 -1e-9 1e-9])
             drawnow
             subplot(1,4,2)
             hold on
            plot(Rs,hguess1*hinitial,'red')
            plot(Rs,hcalc*hinitial,'blue')
            axis([0 100e-6 0e-9 3000e-9])
           title(num2str(runs))
         
             drawnow
            if criteriacount > 0 || runs <6 || indexer < 11
                runs=runs+1;
            end
            if runs > 1000000
                break
            end
      
        end
        
    end
  cla(subplot(1,4,2))
 cla(subplot(1,4,3))
   subplot(1,4,1)
    hold on
   plot(Rs,hcalc*hinitial)
    axis([0 250e-6 0 3500e-9])
    title(num2str(indexer));
    drawnow

    VDW_force = trapz(Rsdless,2*pi*Rsdless.*disjoining_p_VDW);
    DL_force = trapz(Rsdless,2*pi*Rsdless.*disjoining_p_DL);
    
    
    allforce(indexer,1)=force;
    allhvst(indexer,:)=hcalc;
    allpvst(indexer,:)=p1guess;
    allwvst(indexer,:)=wguess;
    allVDW(indexer,:)=VDW_force;
    allDL(indexer,:)=DL_force;
    allxvst(indexer,:)=xxx;
    allruns(indexer,:)=runs;
    time=time+tincredless;
    
    horzcat(alldata, indexer);
    
    qws = num2str(thicknessmax);
  
        qwsstring = strcat(qws,'_', num2str(thicknessdless), '_', num2str(R),'_', num2str(k), '_', num2str(ShiftLength), 'tincre0.5V61h030shift10nm','.xlsx');
        
    hvststring = strcat(qwsstring,'_','allhvst','.csv');
		pvststring = strcat(qwsstring,'_','allpvst','.csv');
		wvststring = strcat(qwsstring,'_','allwvst','.csv');
		forcestring = strcat(qwsstring,'_','allforce','.csv');
		runstring = strcat(qwsstring,'_','allruns','.csv');

	csvwrite(hvststring,allhvst)
	csvwrite(pvststring,allpvst)
	csvwrite(wvststring,allwvst)
	csvwrite(forcestring,allforce)
	csvwrite(runstring,allruns)
	


    indexer=indexer+1;

    if runs > 1000000
        break
    end
     if indexer > 1000
        break
    end
end





end


