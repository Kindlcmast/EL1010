%% Lab 3 Johan Kindlundh, KINDL
clear all
clc
close all

%% Givna konstanter
figureIndex=1;      %För figurer
YYMMDD = 010715;    % ger NAN med 000715
e0=0;               % statiskt fel
ramperror=0.05;     % Ramperror 
Lm = 2;             % Induction
Rm = 21;            % Resistance
b = 1;              % Friction coefficient
Ktau = 38;          % Material constant
Km = 0.5;           % Material constant
n = 1/20;           % Gearing factor 
[J,umax] = lab3robot(YYMMDD);       % J = Moment of inertia, umax = The maximum value of the control signal
%% _________________Assignment 1___________________

%Transfer function G(s)
disp('____Assignment 1____')
s = tf('s');
G = (Ktau*n)/(J*Lm*s^3+(Rm*J+Lm*b)*s^2+(Km*Ktau+b*Rm)*s);
%% %________________Assignment 2__________________
disp('____Assignment 2____')
Kp1 = 8; 
Fp = Kp1;
Gc = (G*Fp)/(1+(G*Fp)); %Closed loop system
disp('Stegsvar:')
stepinfo(Gc)
%Gc=(Ktau*n*Kp)/((Kp*b*Lm*(s^2))+(Kp*J*Lm*(s^3))+(Kp*J*Rm*(s^2))+(Ktau*Km*s)+(Kp*b*Rm*s));
%% ________________Assignment 3__________________
[~,Pm1,Wcg,Wcp1] = margin(G*Fp); % Wcp=Cross over freq, Pm=Phase margin Assignment 3
disp('____Assignment 3____')
disp(['crossover frequency: ' num2str(Wcp1)])
disp(['Phase margin: ' num2str(Pm1)])
disp(['Bandwith of the closed loop system ' num2str(bandwidth(G*Fp/(1+G*Fp)))])
disp(['Peak gain closed loop system: ' num2str(20*log10(getPeakGain(G*Fp/(1+G*Fp))))])
%% _________________Assignment 4___________________  
disp('____Assignment 4____')
%______Plottar_____
figure(figureIndex);figureIndex=figureIndex+1;
margin(G); hold on;
margin(G*Fp); legend('G','GFp'); grid on; hold off;
Kpmax = abs(evalfr(G,Wcg*sqrt(-1)))^-1;
disp(['maximala Kp: ' num2str(Kpmax)])

figure(figureIndex);figureIndex=figureIndex+1;
nyquistplot(G*Fp, G*Kpmax); hold on; title('Nyquist för P-regulator öppet system');grid on; legend('GFp', 'G*Kpmax')
xlim([-1.3,0.2]); ylim([-4,4]);
hold off;
%% _________________Assignment 5___________________

%____________Konstruktion av Lead-Lag____________

%________________Lead___________________
W = Wcp1*4;                           %Önskad skärfrekvens
Kpi = abs(evalfr(G,W*sqrt(-1)))^-1;    %för att få ut fasmarginal vid den skärfrekvensen och för utformning av Lag-controller
[~,Pm2,~,Wcp2] = margin(G*Kpi);        %Kollar hur myket fasen  måste ökas vid skärfrekvensen
lagp = 5;                             %ungefär vad man kommer tappa på lag-delen 
p = Pm1-Pm2+lagp;                     %önskad fasanvancering
%____Beräkning av beta____
x = tand(p);sec = cosd(p)^-1;
beta = 0.5*(4*(x^2)-4*x*sec+2);       % Beta->(0,1) 
%______________________________
Kp = ((1/(beta^0.5))*abs(evalfr(G,W*sqrt(-1))))^-1;     % Kp 
td = (W*sqrt(beta))^(-1);                               % Väljer td så att fasökningen sker vid nya skärfrekvensen
FLead = (td*s+1)/((td*s*beta)+1);                       %lead 
Fl = Kp*FLead;                                          %P och lead
%Kpi = (abs(evalfr(G,W*sqrt(-1))))^-1;                   


%______Plottar_____
figure(figureIndex);figureIndex=figureIndex+1; subplot(2,1,1);
bode(G*Fp,Fl*G,FLead);hold on;grid on
legend('G*Fp','G*Fl','FLead');
[~,Pm3,~,Wcp3] = margin(G*Fl); hold off;
disp('____Assignment 5____')
disp('_______Flead______')
disp(['cross over frequency: ' num2str(Wcp3)])
disp(['Phase margin: ' num2str(Pm3)])

%________________LAG____________________
%1/(lim s->0 (sFsGs))=e1 => lim s->0 Fs = 421.0530
%s->0 Fd = Kp;
%s-> 0 Flag=421.0530/Kp => gamma=kp/421.0530;

gamma = Kp/(((Km*Ktau+b*Rm)/(Ktau*n))/ramperror);
ti = 12.3/W;                        %tummregel ish 1/t1=wc/10
FLag = (ti*s+1)/((ti*s)+gamma);     %lag
F=(Fl*FLag);                        %Fullständiga regulatorn
%______Plottar_____
subplot(2,1,2);
bode(G*F,Kpi*FLag); hold on; grid on;
legend('G*F','Kpi*FLag'); hold off;
[~,Pm4,~,Wcp4] = margin(G*(Kpi*FLag));

disp('_______Flag______')
disp(['cross over frequency: ' num2str(Wcp4)])
disp(['Phase margin: ' num2str(Pm4)])
%_______________________________________

[~,Pm5,~,Wcp5] = margin(G*F);

disp('_______LeadLag______')
disp(['cross over frequency: ' num2str(Wcp5)])
disp(['Phase margin: ' num2str(Pm5)])
disp(['Bandwith of the closed loop system ' num2str(bandwidth(G*F/(1+G*F)))])
disp(['Peak gain closed loop system: ' num2str(20*log10(getPeakGain(G*F/(1+G*F))))])
%_______Ramp response_______
figure(figureIndex);figureIndex=figureIndex+1;
step((G*F/(1+G*F))/s); hold on;
plot(linspace(0,1000),linspace(0.05,1000+0.05),linspace(0,1000),linspace(-0.05,1000-0.05), 'LineWidth',1.5)
xlim([998,1001]);ylim([998,1001]);
legend('ramp','övre gräns', 'undre gräns');title('ramp response LeadLag');hold off;

%_______Figur för utsignal_______
figure(figureIndex);figureIndex=figureIndex+1;
step(F/(1+G*F));hold on; grid on;title('output from LeadLagcontroller at a unit step input'); hold off;
%_______Prestanda_______
disp('Stegsvar:')
stepinfo((G*F/(1+G*F)))
%% _________________Assignment 6___________________
%se dokumentation
%% _________________Assignment 7___________________
%se dokumentation
%% _________________Assignment 8___________________
S1=1/(1+G*Fp);              %känslighetsfunktionen för p-regulatorn
S2=1/(1+G*F);               %känslighetsfunktionen för Leadlag-regulatorn
figure(figureIndex);figureIndex=figureIndex+1;
bodemag(S1, S2);hold on; grid on; title('Bode plot av känslighetsfunktioerna');
legend('S1','S2'); hold off;

%% _________________Assignment 9___________________
DeltaG1 = (s+10)/(40); %DeltaG från uppgiften
DeltaG2 = (s+10)/(4*s+0.04); %DeltaG från uppgiften
T= G*F/(1+G*F);
%om |1/T|>|detaG| kan funktionen anses vara robust, dock kan funktionen anses vara robust trots att kriteriet inte är uppfyllt 
figure(figureIndex);figureIndex=figureIndex+1;
bodemag(DeltaG1,DeltaG2,1/T);hold on; grid on;
legend('deltaG01','deltaG02', '1/T'); hold off;

%% _________________Assignment 10___________________
%Se dokumentation för beräknign
A=[0,n,0;0,-b/J, Ktau/J;0,-Km/Lm,-Rm/Lm];
B=[0;0;Lm^(-1)];
C=[1,0,0];
D=0;
sys=ss(A,B,C,D);
%% _________________Assignment 11___________________
S= [B,A*B,(A^2)*B];                     %Styrbarhetsmatrisen 
Dets=det(S);
%Är determinanten skilld från noll kommer systemet vara styrbart
disp('____Assignment 11____')
if Dets==0
    disp('Systemet är inte styrbart')
else
    disp('Systemet är styrbart')
end
O=[C;C*A;C*A^2];                        %observerbarhetsmatrisen 
DetO=det(O);
%Är determinanten skilld från noll kommer systemet vara observerbart
if DetO==0
    disp('Systemet är inte observerbart')
else
    disp('Systemet är observerbart')
end

%% _________________Assignment 12___________________
%poleplacement
%Lead lag controller: d=0.623 W0=2.45, tp=10.9, och några närmare origo,
%men de försvinner fort pga nollställen
d = 0.6852;                                               % Dämpning 
W0 = 2.5415;                                              % Frekvens 
tp = 10;                                                  % tredje, snabbare polen
p = [-tp,-W0*d+W0*sqrt((d^2)-1),-W0*d-W0*sqrt((d^2)-1)];  % Poler
L = place(A,B,p);                                         % ändrar egenvärdena så att polerna p ges
Acl = A-B*L;                                              % skapar det slutna systemt
syscl = ss(Acl,B,C,D);                                    % skapar det slutna systemt
L0 = (dcgain(syscl))^-1;                                  % väljer L0 så att static gain är noll
syscl = ss(Acl,B*L0,C,D);                                 % skapar det slutna systemt

disp('____Assignment 12____')
disp('________State feedback controller________')
disp(['Bandwith oft the closed loop system ' num2str(bandwidth(syscl))])
disp(['Peak gain closed loop system: ' num2str(20*log10(getPeakGain(syscl)))])
stepinfo(syscl)
figure(figureIndex);figureIndex=figureIndex+1;hold on; grid on; title('step response of state feedback controller')
step(syscl); hold off;

figure(figureIndex);figureIndex=figureIndex+1;
rlocus(syscl,(G*F/(1+G*F))); hold on; grid on;
legend('State feedback controller', 'Lead lag controller'); hold off;

figure(figureIndex);figureIndex=figureIndex+1;hold on; grid on;
bode(syscl,T);legend('State feedback controller', 'Lead lag controller'); hold off;

%% _________________Assignment 13___________________
%se dokumentation
%% _________________Verifiering_____________________
lab3robot(G,Kp1,F,A,B,C,L,L0,YYMMDD) 






