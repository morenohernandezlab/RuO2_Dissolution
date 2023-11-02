clear all
load RuO2_manuscript_video_info_good.mat video_info_good

%% Loot at dimensions
clearvars -except video_info_good
v = 1;

report_fit = [video_info_good(v).time(video_info_good(v).frames),video_info_good(v).fit_all(:,1:7)];
report_dim = [video_info_good(v).time(video_info_good(v).frames),report_fit(:,2)+report_fit(:,3),report_fit(:,4)+report_fit(:,7),report_fit(:,5)+report_fit(:,6),report_fit(:,8)];
offset = 12;

report_dim = report_dim(offset:end,:);
report_dim(:,1) = report_dim(:,1)-report_dim(1,1);
video_info_good(v).offset = offset;
%
hold on
tiledlayout(2,2)
for i = 1:4
    nexttile
plot(report_dim(:,1),report_dim(:,1+i))
end
hold off

%% Fit data
k_rates = linspace(0,3,31);
t_array = [.1,1,2,4]
s_array = [2,2.1,2.2,2.3,2.4,2.5,2.75,3,4]
[ka,ta,sa] = ndgrid(k_rates,t_array,s_array);
trial_values = [ka(:),ka(:),ka(:),ka(:),ta(:),sa(:)];
time_end = report_dim(end,1); %Seconds
M_RuO2 = 52.4; %Molar density (mol/L) of RuO2 material
M_s = 19.9; %Molar density (mol/L) of saturated solution
tol = 4;
opts = odeset('AbsTol',10^-tol,'RelTol',10^-tol);
dim_initial = report_dim(1,2:5);
V_initial = (2*dim_initial(1))*(2*dim_initial(1))*dim_initial(4); %volume(shp3d)



for i = 1:size(trial_values,1) %make par
i
k_rates = trial_values(i,1:4); %Reaction rate for etchin in nm/s
t_rad = trial_values(i,5); %Radiolysis time constant in s, set max to 10 s
s = trial_values(i,6); % Radius of liquid pocket in nm, usually 1.1 to 10
ode_param = [M_RuO2,M_s,dim_initial,k_rates,t_rad,s,V_initial];
[time,dimensions] = ode45(@(t,y) myode3(t,y,ode_param),[0,time_end], [dim_initial],opts);
dimensions = max(dimensions,0);
interp_dimensions = interp1(time,dimensions,report_dim(:,1));
clf
hold on
plot(report_dim(:,1),interp_dimensions(:,1))
plot(report_dim(:,1),report_dim(:,2))

hold off 
drawnow
error_fit(i,:) = [i,sum(((interp_dimensions(:,[1])-report_dim(:,2))).^2,'all')];


end
%
error_no_fail = [];
for i = 1:size(error_fit,1)
if error_fit(i,1) ~= 0
error_no_fail = [error_no_fail;error_fit(i,:)];
end
end
[~,min_indx] = min(error_no_fail(:,2))
best_fit = trial_values(error_no_fail(min_indx,1),:); %[t_a(:),k_L_a(:),k_W_a(:),s_a(:)];
ode_param = [M_RuO2,M_s,dim_initial,best_fit(1:4),best_fit(6),best_fit(6),V_initial];
[time,dimensions] = ode45(@(t,y) myode3(t,y,ode_param),[0,time_end], dim_initial,opts);
interp_dimensions = interp1(time,dimensions,report_dim(:,1));
clf
subplot(2,2,1)
for k = 1:4
    nexttile
hold on
plot(report_dim(:,1),interp_dimensions(:,k))
plot(report_dim(:,1),report_dim(:,k+1))

hold off 
end
drawnow


% Refine fits
for u = 1:10
for j = 1:4
p = .4;
counts = 0;
while p > .0005
k_rates = linspace(1-p,1+p,9);

trial_values = k_rates(:);
error_fit = [];
for i = 1:size(trial_values,1) %make par
k_rates = best_fit(1:4);
k_rates(j) = k_rates(j)*trial_values(i);
t_rad = best_fit(5);
s = best_fit(6);
ode_param = [M_RuO2,M_s,dim_initial,k_rates,t_rad,s,V_initial];
[time,dimensions] = ode45(@(t,y) myode3(t,y,ode_param),[0,time_end], [dim_initial],opts);
interp_dimensions = interp1(time,dimensions,report_dim(:,1));
error_fit(i,:) = [i,sum(((interp_dimensions(:,j)-report_dim(:,j+1))).^2,'all')];
end
error_no_fail = [];
for i = 1:size(error_fit,1)
if error_fit(i,1) ~= 0
error_no_fail = [error_no_fail;error_fit(i,:)];
end
end

[~,min_indx] = min(error_no_fail(:,2))
best_fit(j) = best_fit(j)*trial_values(min_indx);



if min_indx == 5
    p = p*.9; %/2
end
counts = counts +1;
if counts > 5
p = p*.9; %/2
end

end

end
% See best fit
ode_param = [M_RuO2,M_s,dim_initial,best_fit(1:4),best_fit(5),best_fit(6),V_initial];

[time,dimensions] = ode45(@(t,y) myode3(t,y,ode_param),[0,time_end], dim_initial,opts);
%dimensions = max(dimensions,0);
interp_dimensions = [];
interp_dimensions = interp1(time,dimensions,report_dim(:,1));
clf
for k = 1:4
    nexttile
hold on
plot(report_dim(:,1),interp_dimensions(:,k))
plot(report_dim(:,1),report_dim(:,k+1))

hold off 
end
drawnow
end

% save best fit
video_info_good(v).best_fit = best_fit;



%% View frames and define tip type
clearvars -except video_info_good
v = 37
frames = video_info_good(v).frames;

frame_snapshot = [frames(1),frames(end/2),frames(end)];
clf
subplot(1,3,1)
for i = 1:3
    nexttile
imshow([video_info_good(v).Ibw{frame_snapshot(i)}])

end
%video_info_good(v).type = '111'
%video_info_good(v).type = '001'
%video_info_good(v).type2 = 'both'

%% Label Control Nanocyrstals

%video_info_good(v).control = 'yes'

%% Save labeling of dataset

save RuO2_manuscript_video_info_good_controllabelled_v2.mat video_info_good

%% create structures of the wires
clear all
load RuO2_manuscript_video_info_good_controllabelled_v2.mat video_info_good

tips_111 = [];
tips_001 = [];
tips_both = [];
for i = 1:37
if and(isequal(video_info_good(i).type,'111'),isequal(video_info_good(i).control,[]))
tips_111 = [tips_111; [i,video_info_good(i).fit_all(1,1:7),video_info_good(i).best_fit]];

end
if and(isequal(video_info_good(i).type,'001'),isequal(video_info_good(i).control,[]))
tips_001 = [tips_001; [i,video_info_good(i).fit_all(1,1:7),video_info_good(i).best_fit]];

end
if and(isequal(video_info_good(i).type2,'both'),isequal(video_info_good(i).control,[]))
tips_both= [tips_both; [i,video_info_good(i).fit_all(1,1:7),video_info_good(i).best_fit]];

end
end

%% View statistics for 001 tipped wires

tip_001_110_rate = tips_001(:,9);
tip_001_001_rate = tips_001(:,12);
tip_001_001_over_111_rate = tip_001_001_rate./tip_001_110_rate;
clf
hold on
plot(tip_001_110_rate,tip_001_001_rate,'o')
plot([0,5],[0,5]*1)

hold off
xlabel('(110) etching rate (nm s^{-1})')
ylabel('(001) etching rate (nm s^{-1})')

%% View statistics for 111 tipped wires
tip_111_110_rate = tips_111(:,9);
tip_111_111_rate_1 = tips_111(:,10);
tip_111_111_rate_2 = tips_111(:,11);

tip_111_111_over_111_rate = .5*(tip_111_111_rate_1+tip_111_111_rate_2)./tip_111_110_rate;

clf
hold on
plot(tip_111_110_rate,tip_111_111_rate_1,'o')
%plot(tip_111_111_rate_1,tip_111_111_rate_2,'o')
plot([0,10],[0,10]*1)
hold off
xlabel('(110) etching rate (nm s^{-1})')
ylabel('(111) etching rate (nm s^{-1})')


%% THis section of the code is used to determine the relative rates with the dimension vs. dimension method
load RuO2_manuscript_video_info_good_controllabelled.mat
clearvars -except video_info_good
clf
for v = 7

report_fit = [video_info_good(v).time(video_info_good(v).frames),video_info_good(v).fit_all(:,1:7)];
report_dim = [video_info_good(v).time(video_info_good(v).frames),report_fit(:,2)+report_fit(:,3),report_fit(:,4)+report_fit(:,7),report_fit(:,5)+report_fit(:,6),report_fit(:,8)];
%report_dim(end-4:end,:) = [];
%report_dim(1,:) = [];

subplot(3,3,1)
hold on
x1 = report_dim(1,2)-report_dim(:,2);
y1 = report_dim(1,3)-report_dim(:,3);
plot(x1,y1,'o')
c1 = fitlm(x1,y1);
video_info_good(v).fit_111_01 = [c1.Coefficients.Estimate(1),c1.Coefficients.Estimate(2),c1.Rsquared.Ordinary];

plot(x1,c1.Coefficients.Estimate(2)*x1+c1.Coefficients.Estimate(1))
hold off
subplot(3,3,2)
plot(report_dim(:,1),report_dim(:,2))
subplot(3,3,3)
plot(report_dim(:,1),report_dim(:,3))
subplot(3,3,4)
hold on
y2 = report_dim(1,4)-report_dim(:,4);

plot(x1,y2,'o')
c2 = fitlm(x1,y2);
video_info_good(v).fit_111_02 = [c2.Coefficients.Estimate(1),c2.Coefficients.Estimate(2),c2.Rsquared.Ordinary];
plot(x1,c2.Coefficients.Estimate(2)*x1+c2.Coefficients.Estimate(1))
hold off
subplot(3,3,5)
plot(report_dim(:,1),report_dim(:,2))
subplot(3,3,6)
plot(report_dim(:,1),report_dim(:,4))



subplot(3,3,7)
hold on
y3 = report_dim(1,5)-report_dim(:,5);

plot(x1,y3,'o')
c3 = fitlm(x1,y3);
video_info_good(v).fit_001 = [c3.Coefficients.Estimate(1),c3.Coefficients.Estimate(2),c3.Rsquared.Ordinary];

plot(x1,c3.Coefficients.Estimate(2)*x1+c3.Coefficients.Estimate(1))
hold off
subplot(3,3,8)
plot(report_dim(:,1),report_dim(:,2))
subplot(3,3,9)
plot(report_dim(:,1),report_dim(:,5))

%video_info_good(v).fit_dim_vs_dim = [c1; c2; c3];
end


%%

clear all
load RuO2_manuscript_video_info_good_controllabelled_v2.mat video_info_good

%% Pull up data to be plotted
clear all
load RuO2_manuscript_video_info_good_controllabelled_v2.mat

%% select video
v = 12
best_fit = video_info_good(v).best_fit
offset = 1%video_info_good(v).offset;
report_fit = [video_info_good(v).time(video_info_good(v).frames),video_info_good(v).fit_all(:,1:7)];
report_dim = [video_info_good(v).time(video_info_good(v).frames),report_fit(:,2)+report_fit(:,3),report_fit(:,4)+report_fit(:,7),report_fit(:,5)+report_fit(:,6),report_fit(:,8)];
dim_initial = report_dim(offset,2:5);
M_RuO2 = 52.4; %Molar density (mol/L) of RuO2 material
M_s = 19.9; %Molar density (mol/L) of saturated solution
tol = 6;
opts = odeset('AbsTol',10^-tol,'RelTol',10^-tol);
k_rates = best_fit(1:4)
t_rad = best_fit(5)
s = best_fit(6);
time_end = report_dim(end,1);
V_initial = (2*dim_initial(1))*(2*dim_initial(1))*dim_initial(4); %volume(shp3d)

ode_param = [M_RuO2,M_s,dim_initial,k_rates,t_rad,s,V_initial];
[time,dimensions] = ode45(@(t,y) myode3(t,y,ode_param),[0,time_end], [dim_initial],opts);
clf

dim_model = [time,dimensions];
dim_plot = 5;
hold on
plot(report_dim(:,1),report_dim(:,dim_plot),'o')
plot(dim_model(:,1),dim_model(:,dim_plot))
%ylim([800,1200])
hold off
%%
function dydt = myode3(t,y,ode_param)


%ode_param = [M_RuO2,M_s,dim_initial,k_rates,t_rad,s];
%ode_param = [M_RuO2,M_s,dim_initial,k_rates,t_rad,s,V_initial];
M_RuO2 = ode_param(1); %Molar density (mol/L) of RuO2 material
M_s = ode_param(2); %Molar density (mol/L) of saturated solution
dim_initial = ode_param(3:6);
k1 = ode_param(7);
k2 = ode_param(8);
k3 = ode_param(9);
k4 = ode_param(10);
t_rad = ode_param(11); %Radiolysis time constant in s
s = ode_param(12); % Radius of liquid pocket in nm, usually 1.1 to 10
V_initial = ode_param(13);
volume_shape = (2*y(1))*(2*y(1))*y(4);
v_ratio =  volume_shape./V_initial;
v_ratio = min(v_ratio,1);
%v_change = max(v_change,0)
V_t = V_initial*v_ratio; %volume(shp3d);

f = (V_initial-V_t)./(s*V_initial-V_t);
%(1-((M_RuO2/M_s)*f))

dydt(1) = -k1*(1-((M_RuO2/M_s)*f))*(1-exp(-t./t_rad)); 
dydt(2) = -k2*(1-((M_RuO2/M_s)*f))*(1-exp(-t./t_rad)); 
dydt(3) = -k3*(1-((M_RuO2/M_s)*f))*(1-exp(-t./t_rad)); 
dydt(4) = -k4*(1-((M_RuO2/M_s)*f))*(1-exp(-t./t_rad)); 


dydt = dydt';
end




