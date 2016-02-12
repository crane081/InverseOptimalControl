%% Estimating individual's sensormotor speed in a simulated driving task
%  Crane Huang, 02/11/2016.

%%---------------------- Experiment ----------------------%%
% Subjects are instructed to push the joystick to the maximal forward
% position as soon as they perceive the car starts moving.

%%----------------------- Variables ----------------------%%
% rt: reaction time from car motion onset to start pushing the joystick
% mt: time used to push the joystick from resting to target position
% v:  car speed in each trial (before subject started to push the joystick)
% action: joystick action in each trial

%%--------------------- Sensory speed --------------------%%
% Car dynamic:    dXt=vdt
% Observation:    dYt=gamma(Xt-Yt)dt
% Reaction time:  RT=argmin(t,Yt>=thd)
% where:
% Xt: car position, v: car speed, gamma: sensory speed

%%----------------------- Motor speed --------------------%%
% Joystick:       dCt=beta(Ut-Ct)dt
% Movement time:  MT=argmin(t,Ct>=Ut)
% where:
% Ut: target joystick position, Ct: joystick position, beta: motor speed


clear
close all
load ExampleData %load example subject's data

%% ----------------Sensory speed----------------%%
%  Car dynamic:    dXt=vdt
%  Observation:    dYt=gamma(Xt-Yt)dt
%  Reaction time:  RT=argmin(t,Yt>=thd)
%  Xt: car position, v: car speed, gamma: sensory speed

Tcm = .9347; %convert pix to cm
dt= .0167;   
T = 5;
s = ceil(T/dt);

v = d1.v*Tcm; %speed in cm/second
rt= d1.rt;
N = length(rt); %number of trials

%search range
s_gamma = .05:.05:2;
s_thd   = .001:.01:.1;

h = waitbar(0,'Please wait...');

for i = 1:length(s_gamma)
    for j = 1:length(s_thd)
        for tr = 1:N %go through each trial
            x = [];
            y = [];
            id = [];
            x(1) = 0;
            y(1) = 0;
            
            for t = 2:s %go through each time point
                dx = v(tr)*dt; 
                x(t) = x(t-1)+dx; %car position  at time t
                dy = s_gamma(i)*(x(t)-y(t-1))*dt;
                y(t) = y(t-1)+dy; %observed car position at time t
            end
            
            id = find(y>=s_thd(j),1); %the first time observation reaches the threshold
            
            if any(id)
                rt_est0(tr) = id*dt; %estimated reaction time
                rt_err(tr) = rt_est0(tr)-rt(tr); %estimation error
                rt_est(i,j,tr) = rt_est0(tr);
            else
                rt_err(tr) = 10;
            end
        end
        err(i,j) = sum(rt_err.^2); %sum error of i,j across all trials
    end 
    waitbar(i/length(s_gamma),h)
end

close(h) %close waitbar
%% find the best gamma and threshold
[minE,ind]=min(err(:));
[id_gamma,id_thd]=ind2sub(size(err),ind);

best_gamma = s_gamma(id_gamma); %find best gamma
best_thd = s_thd(id_thd); %find best thd

rt_e=squeeze(rt_est(id_gamma,id_thd,:)); %best RT estimation

%% Plot data and model prediction
figure;
scatter(v,rt,100,'g','filled');
hold on;
scatter(v,rt_e(1:length(v)),200,'filled','b');
xlim([0 .3]);
set(gca,'fontsize',20);
legend('data','model');
ylabel('reaction time (s) ');
xlabel('car speed (cm/s) ');
title('Reaction Time (all trials)');

%% define interested speed range 
vv = .01:.01:.1; 

m_rt = [];
p_rt = [];
std_rt = [];
for i = 1:length(vv)
    idd = find(abs(v-vv(i))<1e-10);
    m_rt(i) = mean(rt(idd));
    std_rt(i) = std(rt(idd));
    p_rt(i) = mean(rt_e(idd));
end

figure;
scatter(vv,p_rt,200,'filled','b','linewidth',2);
hold on
h1 = shadedErrorBar(vv,m_rt,1.66*std_rt,{'g','linewidth',4},1);
legend('model','data');
xlim([0 .1]);
set(gca,'fontsize',20);
xlabel('car speed (cm/s) ');
ylabel('reaction time (s) ');
title('Reaction time (interested range) ');

%% save sensory speed result
gamma_est = best_gamma;
thd_est = best_thd;
rt_model = p_rt;
rt_data = m_rt;

%% ----------------Motor speed----------------%%
% Joystick:       dCt=beta(Ut-Ct)dt
% Movement time:  MT=argmin(t,Ct>=Ut)
% Ut: target joystick position, Ct: joystick position, beta: motor speed

id = find(d1.mt>.1 & d1.mt<1); %define valid range
mt = d1.mt(id);
N = length(mt); %number of trials
Ut = 20; % target position

ctt = [];
ss=round(max(mt)/dt); %time window

%compile all trials
for j = 1:N    
    k = id(j);
    ct = d1.action{k};
    
    idd = find(abs(ct-Ut)<=.1,1);  %the first time reached the target position
    idt = find(ct>=0);    
    ct = ct(idt(1):idd); %valid action range
  
    if (idd-idt(1)<=ss) 
        if idd-idt(1)<ss
            ct(end+1:ss) = Ut;
        else
            ct = ct(1:ss);
        end
        ctt = [ctt;ct]; 
    end
end

ct_m = mean(ctt);

d_ct = ct_m(2:end)-ct_m(1:end-1);
df_uc = Ut-ct_m';
x1 = df_uc(1:end-1)*dt;
y1 = d_ct;

clear ct
beta = regress(y1',x1)*2;
ct(1)= 0;
id = [];
for t = 2:s
    d_ct = beta*(Ut-ct(t-1))*dt;
    ct(t) = ct(t-1)+d_ct;
end
cu=abs(ct-ut);
id = find(cu<= 1,1);
mt_e = (id-1)*dt;
    
%% save motor speed result
beta_est = beta;
mt_model = mt_e;
mt_data = mean(mt);
display(sprintf('Mean MT: %.2f s; Model MT: %.2f s',mean(mt),mt_model));
