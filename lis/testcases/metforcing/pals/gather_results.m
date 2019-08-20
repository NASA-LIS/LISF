clear all
close all
clc


% direct us to the lis output directory
lis_out_dir = '/discover/nobackup/gnearing/pals_testcase/testcase/output/';

% if you want, you can pre-deignate the number of timesteps
ntimes = 48*366*7.1; 
data = zeros(ntimes,8); 

t=0;
dt = 0;
for y = 2000:2006
  for m = 1:12
    ms = num2str(m);
    if m < 10
      ms = strcat('0',ms);
    end
    for d = 1:31
      useday = 1;
      ds = num2str(d);
      if d < 10
        ds = strcat('0',ds);
      end
      for h = 0:23
        hs = num2str(h);
        if h < 10
          hs = strcat('0',hs);
        end
        for mn = 0:30:30
          if mn == 0
            mns = '00';
          else
            mns = '30';
          end
  
          file = strcat(lis_out_dir,'SURFACEMODEL/',num2str(y),ms,...
                  '/LIS_HIST_',num2str(y),ms,ds,hs,mns,'.d01.nc');

          try 
            Qle = ncread(file,'Qle_tavg');
            Qle = Qle(1);
            Qh = ncread(file,'Qh_tavg');
            Qh = Qh(1);
            Qg = ncread(file,'Qg_tavg');
            Qg = Qg(1);
         
            t = t+1;
            data(t,:) = [y,m,d,h,mn,Qle,Qh,Qg];
      
            fprintf('finsihed year %d, month %d, day %d, hour %d \n',y,m,d,h);
          catch
            useday = 0;
            fprintf('******* year %d, month %d, day %d, hour %d \n',y,m,d,h);
          end % try
        end % minute
      end % hour
    end % day
  end % month
end % year

% clear unused storage
data(data(:,1) == 0,:) = [];

% save results
save('./output/lis_output.txt','data','-ascii');

% gather PALS results
file = './output/Blodgett_PALS_Obs.nc';
Qle = squeeze(ncread(file,'Qle'));
Qh = squeeze(ncread(file,'Qh'));
Qg = squeeze(ncread(file,'Qg'));

% plot results
figure(1); close(1); figure(1);
plot(data(:,6:end));
xlabel('Timestep','fontsize',16);
ylabel('Surface Fluxes [W/m2]','fontsize',16);
legend('Q_l_e','Q_h','Q_g');
title('LIS Results');

% plot comparison
figure(2); close (2); figure(2);
subplot(1,3,1)
title('Qle','fontsize',18);
plot(Qle(1:length(data)),data(:,6),'.');
xlabel('Observed','fontsize',16);
ylabel('Modeled','fontsize',16);

subplot(1,3,2)
title('Qh','fontsize',18);
plot(Qh(1:length(data)),data(:,7),'.');
xlabel('Observed','fontsize',16);
ylabel('Modeled','fontsize',16);

subplot(1,3,3)
title('Qg','fontsize',18);
plot(Qg(1:length(data)),data(:,8),'.');
xlabel('Observed','fontsize',16);
ylabel('Modeled','fontsize',16);

% plot results
figure(3); close(3); figure(3);
plot(Qle(1:length(data)),'b'); hold on;
plot(Qh(1:length(data)),'g'); hold on;
plot(Qg(1:length(data)),'r'); hold on;
xlabel('Timestep','fontsize',16);
ylabel('Surface Fluxes [W/m2]','fontsize',16);
legend('Q_l_e','Q_h','Q_g');
title('LIS Results');






