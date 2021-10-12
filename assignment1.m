%% SETUP AND SETTINGS
cleansheet
% Na current
f_ph_Na=0.5;
Na_trend = linspace(1,f_ph_Na, 15);
% Na_trend = Na_trend(1:10);

% CaL current
f_ph_CaL=0.5;
CaL_trend = linspace(1,f_ph_CaL, 15);
% CaL_trend = CaL_trend(1:10);

% NaCa current
f_ph_NaCa=0.1;
NaCa_trend = linspace(1,f_ph_NaCa, 15);
% NaCa_trend = NaCa_trend(1:10);

settings.NumStim = 75;
settings.storeLast = 5;

%% SIMULATIONS
max_minutes = 15;
Vr = zeros(max_minutes,1);
Vpp = zeros(max_minutes,1);
dVdTmax = zeros(max_minutes,1);
APD90 = zeros(max_minutes,1);
ffreq = zeros(max_minutes,1);
Vrecord = cell(max_minutes,1);
Trecord = cell(max_minutes,1);

for m = 1:max_minutes
    [StateVars,Ti]=TT_endo_Ko_main(settings,Na_trend(m),CaL_trend(m),NaCa_trend(m));
    V = StateVars(:,10);
    Vrecord{m} = V;
    Trecord{m} = Ti;
    Vr(m) = V(end);
    Vpp(m) = max(V)-min(V);
    dVdTmax(m) = max(diff(V)./[0 Ti]);
    thr_90_rep = max(V)-0.9*(max(V)-V(end));
    rep_thr_crossed = Ti(abs(V-thr_90_rep)<2 & diff([0; V])<-0.1);
    APD90(m) = rep_thr_crossed(1);
    peaks_idx = [];
    for i=2:length(V)-1
        if V(i)>V(i-1) && V(i) > V(i+1)
           peaks_idx = cat(1,peaks_idx,i); 
        end
    end
    peaks_idx = peaks_idx(1:2:end);
    ffreq(m) = 1/(Ti(peaks_idx(2))-Ti(peaks_idx(1)))*1000;
end

%% PLOTTING
figure('Name','AP minute-wise');
minute = 1:max_minutes;
for m=minute
    plot(Trecord{m},Vrecord{m});
    hold on;
end
xlabel("Time [msec]",'interpreter','latex');
ylabel("$V_m$ [mV]",'interpreter','latex');
hold off; grid on; title('\textbf{APs over time}','interpreter','latex');
legend('min1','min2','min3','min4','min5','min6','min7','min8','min9',...
    'min10','min11','min12','min13','min14','min15','interpreter','latex','Location','eastoutside');

figure('Name','Biomarkers');
subplot(221);
plot(minute,Vr,'*-','LineWidth',2); grid on; 
title("\textbf{Rest Potential}",'interpreter','latex')
xlabel("Minute after occlusion [min]",'interpreter','latex'); 
ylabel("$V_{rest}$ [mV]",'interpreter','latex'); 
xlim([0 10]);

subplot(222);
plot(minute,Vpp,'*-','LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]); grid on; 
title("\textbf{Peak-to-Peak voltage}",'interpreter','latex')
xlabel("Minute after occlusion [min]",'interpreter','latex');
ylabel("$V_{pp}$ [mV]",'interpreter','latex'); 
xlim([0 10]);

plot(minute,dVdTmax,'*-','LineWidth',2,'Color',[0.4660, 0.6740, 0.1880]	); grid on; 
title("\textbf{Maximum depolarization slope}",'interpreter','latex')
xlabel("Minute after occlusion [min]",'interpreter','latex'); 
ylabel("$\frac{dV}{dt}_{max}$ [mV/ms]",'interpreter','latex'); 
xlim([0 10]);
subplot(224);
plot(minute,APD90,'*-','LineWidth',2,'Color',[0.6350, 0.0780, 0.1840]); grid on; 
title("\textbf{APD90}",'interpreter','latex')
xlabel("Minute after occlusion [min]",'interpreter','latex'); 
ylabel("$APD_{90}$ [msec]",'interpreter','latex'); 
xlim([0 10]);

% % Answer to this question giving your opinion about: Acidosis in 
% % Cardiomyocytes could have a pro-arrhythmic effect?
