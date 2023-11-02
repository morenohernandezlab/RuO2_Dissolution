clear all
clc
tol_used = 10^-8; %Model tolerance
time_run = 10^3;
beam_dose = 2000;
%Composition: MgCl2*6H2O and FeCl3*6H2O
%Density of MgCl2*6H2O: 1.57 g/cm3
%Density of FeCl3*6H2O: 1.82 g/cm3

f_w = 1; %(12)/(6+1+2+6+1+3); %Water fraction
f_cl = 0; %(5)/(6+1+2+6+1+3); %Chloride fraction

Fe_range = linspace(0,1,101);
Cl_range = linspace(0,1,101);


[Fe_map, Cl_map] = meshgrid(Fe_range, Cl_range);
trial_map = [Fe_map(:),Cl_map(:)];

% define species names

names = ["electrons";
    "H_ions";
    "OH_ions";
    "H2O2";
    "HO2_ions";
    "H_mono";
    "OH";
    "O_ions";
    "HO2";
    "O2_ions";
    "O2";
    "H2";
    "O3_ions"
    "O3";
    "HO3";
    "H2O";
    "Cl_ions";
    "ClOH_ions";
    "Cl_mono";
    "Cl2_ions";
    "Cl3_ions";
    "Cl2";
    "Fe2_ions";
    "Fe3_ions";
    "Fe4_ions";
    "Fe5";
    "Fe6_ions";
    "Fe3_OH_ions";
    "Fe3_OH2_ions";
    "Fe3_dimer_ions";
    "Fe2_OH_ions";
    "Fe3_HO2_ions";
    "Fe3_H2O3_ions";
    "O2_2m_ions";
    "Fe2_Cl_ions";
    "Fe3_Cl_ions";
    "Fe3_Cl2_ions";
    ];

tic
parfor i = 1:size(trial_map,1)
    i
salt = 100;
pH_initial = 1;
H2O_initial = 55.6; %(Fe_salt_molarity*6+Mg_salt_molarity*6); %M, assuming salts are FeCl3*6H2O and MgCl2*6H2O
Cl_initial = trial_map(i,2); %salt*0*10^-3+100*10^-3; %(Fe_salt_molarity*3+Mg_salt_molarity*2); %M
Fe2_initial = 0;
Fe3_initial = trial_map(i,1); %0.*(Fe_salt_molarity*1); %M
Fe4_initial = 0;
Fe5_initial = 0;
Fe6_initial = 0;
Fe_rates =[2.7*10^8 2.7*10^6, 2.0*10^10]; %OH, H_mono, electrons


initialConcentration = [0 10^-pH_initial 10^-(14-pH_initial) 0 0 0 0 0 0 0 0 0 0 0 0 H2O_initial Cl_initial 0 0 0 0 0 Fe2_initial Fe3_initial Fe4_initial Fe5_initial Fe6_initial 0 0 0 0 0 0 0 0 0 0]*10^6;
timeSpan = [10^-8,time_run];
%
%Define constants
numberOfSpecies = length(initialConcentration);

%Molar to micromolar conversion factor%
molarToMicromolar = 10^6;


%Set of Dose Rate in Gy/s (A)%
SP = 2.798; %2.798 MeVcm2/g for 200 keV, 2.360 MeVcm2/g for 300 keV
conversion_factor = SP*10^16*1.60218*10^-13*10^3;
doseRate = beam_dose*conversion_factor;%Conversion of dose rate (e- A-2 s-1) for 300 KeV electrons to Gy/s;



% G-Values (molecules/100 eV)%
% Default G-Values taken from Hill and Smith at 300KeV [ref] %
 gValues = zeros(numberOfSpecies,1);
 gValues(1)  = 3.47 / 100;  % Hydrated Electrons
 gValues(2)  = 4.42 / 100;  % H+
 gValues(3)  = 0.95 / 100;  % OH-
 gValues(4)  = 0.47 / 100;  % H2O2
 gValues(6)  = 1.00 / 100;  % H
 gValues(7)  = 3.63 / 100;  % OH
 gValues(12) = 0.17 / 100;  % H2
 %gValues(12) = conditions_test(i,1) / 100;  % H2 lowest yield reported by LaVerne 2005
 gValues(9) = 1/3*(gValues(2)-gValues(3)-2*gValues(4)+gValues(6)-gValues(7)+2*gValues(12));% HO2
 gValues(16) = -gValues(7)-2*gValues(4)-2*gValues(9) - gValues(3); %Hydrogen balance
 g_O = - gValues(3)-2*gValues(4)-gValues(7)-2*gValues(9) ; %Oxygen balance
 g_H = -0.5*(gValues(2)+gValues(3)+2*gValues(4)+gValues(6)+gValues(7)+2*gValues(12)+gValues(9)); %Hydrogen balance balance
 difference = g_O-g_H;
 
 % Correct for water fraction
gValues = gValues*f_w;

%Direct chloride ionization

gValues(19) = 4.42 / 100; %Cl* generation;
gValues(19) = gValues(19)*f_cl;
gValues(1)  = gValues(1)+gValues(19);
gValues(17) = gValues(17)-gValues(19);

%avogadro's number%
avogadro = 6.022*10^23;

%convert from molecules/(100 eV) to molecules/eV
% gValues = gValues./100;

%conversion from cubic meters to liters%
m3_L = 1000;

%%calculated in micromolar/s
doseRate = doseRate / (1.6e-22);
generationDueToGValues = molarToMicromolar * doseRate.*gValues./(avogadro*m3_L);
 % Define rate constants %
%equilibria%
K = zeros(5,1);

%H2O <=> H+ + OH-&
K(1) = 10^-13.999; %Molar^2

%H2O2 <=> H+ + HO2-%
K(2) = 10^-11.65;%Molar

%OH <=> H+ + O-%
K(3) = 10^-11.9;%Molar

%HO2 <=> H+ + O2-%
K(4) = 10^-4.57;%Molar

%H <=> H+ + e-%
K(5) = 10^-9.77;%Molar

%Rate constants%
%%1/(Ms) unless specified%%



% Note, some are out of numerical order so they are defined when they are
% used in other rate constants. To convert numbering scheme from what is
% presented here to what it presented in the paper [ref], add 6 to the
% index. (i.e., Reaction 25's rate constant has an index of 19)
k = zeros(73,1);
k(1) = 1.4*10^11;
k(2) = k(1)*K(1)*molarToMicromolar^2;%M/s
k(4) = 5*10^10;
k(3) = k(4)*K(2)*molarToMicromolar;%/s
k(5) = 1.3*10^10;
k(6) = k(5)*(K(1)/K(2))*molarToMicromolar;%%/s
k(7) = 1.9*10;
k(8) = 2.2*10^7;
k(10) = 2.3*10^10;
k(9) = k(10)*K(5)*molarToMicromolar;%/s
k(11) = 1.3*10^10;
k(12) = k(11)*(K(1)/K(3))*molarToMicromolar;%/s
k(14) = 10^11;
k(13) = k(14)*K(3)*molarToMicromolar;%/s
k(16) = 5*10^10;
k(15) = k(16)*K(4)*molarToMicromolar;%/s
k(17) = 5*10^10;
k(18) = k(17)*(K(1)/K(4))*molarToMicromolar;%/s%
k(19) = 3*10^10;
k(20) = 1.1*10^10;
k(21) = 1.3*10^10 ;
k(22) = 2*10^10;
k(23) = 1.9*10^10;
k(24) = 5.5*10^9; 
k(25) = 2.5*10^10 ;
k(26) = 3.5*10^9;
k(27) = 2.2*10^10 ;
k(28) = 1.6*10^10; 
k(29) = 3.6*10^10;
k(30) = 1.1*10;
k(31) = 10^10;
k(32) = 9*10^7;
k(33) = 10^10;
k(34) = 7.8*10^9;
k(35) = 7.0*10^9;
k(36) = 9*10^7;
k(37) = 2.1*10^10;
k(38) = 1.8*10^10;
k(39) = 1.8*10^10;
k(40) = 3.8*10^10;
k(41) = 3.6*10^9;
k(42) = 6*10^9; 
k(43) = 8.2*10^9;
k(44) = 4.3*10^7;
k(45) = 2.7*10^7;
k(46) = 2.5*10^10;
k(47) = 7.5*10^9;
k(48) = 2.6*10^9;
k(49) = 6*10^9;
k(50) = 1.1*10^8;
k(51) = 8*10^7;
k(52) = 7*10^5;
k(53) = 6*10^9;
k(54) = 5*10^-1;
k(55) = 5*10^-1;
k(56) = 6*10^9;
k(57) = 5*10^8;
k(58) = 10^2;
k(59) = 6*10^8;
k(60) = 1.3*10^-1;
k(61) = 1.3*10^-1;
k(62) = 10^4 ;
k(63) = 1.5*10^9;
k(64) = 10^9 ;
k(65) = 3.6*10^9;
k(66) = 8*10^7;
k(67) = 5*10^8;
k(68) = 4*10^8;
k(69) = 7*10^8;
k(70) = 5*10^9;
k(71) = 3.3*10^3*molarToMicromolar;%/s
k(72) = 9*10^10;
k(73) = 1.1*10^5*molarToMicromolar;%/s
k(74) = 4.3*10^9;
k(75) = 6.1*10^9*molarToMicromolar;%/s
k(76) = 8.5*10^9;
k(77) = 2.1*10^10;
k(78) = 2.3*10^1*molarToMicromolar; %/s
k(79) = 1.8*10^10;
k(80) = 6*10^4*molarToMicromolar; %/s
k(81) = 2*10^9;
k(82) = 6.3*10^8;
k(83) = 1.0*10^4;
k(84) = 5.0*10^4*molarToMicromolar; %/s
k(85) = 8.8*10^7;
k(86) = 1.0*10^10;
k(87) = 1.0*10^10;
k(88) = 3.0*10^10;
k(89) = 1.0*10^10;
k(90) = 8.0*10^9;
k(91) = 1.0*10^10;
k(92) = 1.0*10^9;
k(93) = 1.0*10^9;
k(94) = Fe_rates(1);
k(95) = Fe_rates(1);
k(96) = Fe_rates(1);
k(97) = Fe_rates(1);
k(98) = Fe_rates(2);
k(99) = Fe_rates(2);
k(100) = Fe_rates(2);
k(101) = Fe_rates(2);
k(102) = Fe_rates(3);
k(103) = Fe_rates(3);
k(104) = Fe_rates(3);
k(105) = Fe_rates(3);
k(106) = 1*10^6; %
k(107) = 0.9*10^7; %
k(108) = 0.6*10^6; %
k(109) = 1.0*10^3; %
k(110) = 1.0*10^4; %
k(111) = 2.34*10^7*molarToMicromolar;%/s
k(112) = 1.0*10^10;
k(113) = 4.68*10^3*molarToMicromolar;%/s
k(114) = 1.0*10^10*(1/molarToMicromolar); %1/(M^2 *s)
k(115) = 1.12*10^7;
k(116) = 1.0*10^10*(1/molarToMicromolar); %1/(M^2 *s)
k(117) = 1.9*molarToMicromolar;%/s
k(118) = 1.0*10^10;
k(119) = 55;
k(120) = 5.9*10^6;
k(121) = 3.1*10^7;
k(122) = 1.0*10^10;
k(123) = 2.0*10^6;
k(124) = 1.0*10^10;
k(125) = 2.3*10^-3*molarToMicromolar;%/s
k(126) = 2.3*10^-3*molarToMicromolar;%/s
k(127) = 2.7*10^8;
k(128) = 1.2*10^6;
k(129) = 1.0*10^7;
k(130) = 1.0*10^10;
k(131) = 2.0*10^4;
k(132) = 5.0*10^7;
k(133) = 2.88*10^10;
k(134) = 1.0*10^10*molarToMicromolar;%/s
k(135) = 6.61*10^10;
k(136) = 1.0*10^10*molarToMicromolar;%/s
k(137) = 1.05*10^11*(1/molarToMicromolar); %1/(M^2 *s)
k(138) = 1.0*10^10*molarToMicromolar;%/s
k(139) = 55;
k(140) = 5.9*10^9;
k(141) = 5.0*10^6;
k(142) = 1.2*10^6;
k(143) = 1.0*10^7;
k(144) = 2.0*10^4;
k(145) = 2.0*10^4;
k(146) = 5.0*10^7;
k(147) = 5.0*10^7;
k(148) = Fe_rates(1);
k(149) = Fe_rates(1);
k(150) = Fe_rates(1);
k(151) = 6.0*10^9;
k(152) = 6.0*10^9;
k(153) = Fe_rates(3);
k(154) = Fe_rates(3);


conservation_matrix = [];

conservation_matrix(1,:) = [-1, 0, 0, 0, 0]; %electrons = reactants(1); 
conservation_matrix(2,:) = [1, 1, 0, 0, 0]; %H_ions    = reactants(2);
conservation_matrix(3,:) = [-1, 1, 1, 0, 0]; %OH_ions   = reactants(3);
conservation_matrix(4,:) = [0, 2, 2, 0, 0]; %H2O2      = reactants(4);
conservation_matrix(5,:) = [-1, 1, 2, 0, 0]; %HO2_ions  = reactants(5);
conservation_matrix(6,:) = [0, 1, 0, 0, 0]; %H_mono    = reactants(6);
conservation_matrix(7,:) = [0, 1, 1, 0, 0]; %OH        = reactants(7);
conservation_matrix(8,:) = [-1, 0, 1, 0, 0]; %O_ions    = reactants(8);
conservation_matrix(9,:) = [0, 1, 2, 0, 0]; %HO2       = reactants(9);
conservation_matrix(10,:) = [-1, 0, 2, 0, 0]; %O2_ions   = reactants(10);
conservation_matrix(11,:) = [0, 0, 2, 0, 0]; %O2        = reactants(11);
conservation_matrix(12,:) = [0, 2, 0, 0, 0]; %H2        = reactants(12);
conservation_matrix(13,:) = [-1, 0, 3, 0, 0]; %O3_ions   = reactants(13);
conservation_matrix(14,:) = [0, 0, 3, 0, 0]; %O3        = reactants(14);
conservation_matrix(15,:) = [0, 1, 3, 0, 0]; %HO3       = reactants(15);
conservation_matrix(16,:) = [0, 2, 1, 0, 0]; %H2O       = reactants(16);
conservation_matrix(17,:) = [-1, 0, 0, 1, 0]; %Cl_ions   = reactants(17);
conservation_matrix(18,:) = [-1, 1, 1, 1, 0]; %ClOH_ions = reactants(18);
conservation_matrix(19,:) = [0, 0, 0, 1, 0]; %Cl_mono   = reactants(19);
conservation_matrix(20,:) = [-1, 0, 0, 2, 0]; %Cl2_ions  = reactants(20);
conservation_matrix(21,:) = [-1, 0, 0, 3, 0]; %Cl3_ions  = reactants(21);
conservation_matrix(22,:) = [0, 0, 0, 2, 0]; %Cl2       = reactants(22);
conservation_matrix(23,:) = [2, 0, 0, 0, 1]; %Fe2_ions  = reactants(23);
conservation_matrix(24,:) = [3, 0, 0, 0, 1]; %Fe3_ions  = reactants(24);
conservation_matrix(25,:) = [-1, 3, 4, 0, 1]; %Fe4_ions  = reactants(25);
conservation_matrix(26,:) = [0, 3, 4, 0, 1]; %Fe5       = reactants(26);
conservation_matrix(27,:) = [1, 3, 4, 0, 1]; %Fe6_ions  = reactants(27);
conservation_matrix(28,:) = [2, 1, 1, 0, 1]; %Fe3_OH_ions =  reactants(28);
conservation_matrix(29,:) = [1, 2, 2, 0, 1]; %Fe3_OH2_ions =  reactants(29);
conservation_matrix(30,:) = [4, 2, 2, 0, 2]; %Fe3_dimer_ions = reactants(30);
conservation_matrix(31,:) = [1, 1, 1, 0, 1]; %Fe2_OH_ions    = reactants(31);
conservation_matrix(32,:) = [2, 1, 2, 0, 1]; %Fe3_HO2_ions = reactants(32);
conservation_matrix(33,:) = [1, 2, 3, 0, 1]; %Fe3_H2O3_ions = reactants(33);
conservation_matrix(34,:) = [-2, 0, 2, 0, 0]; %O2_2m_ions = reactants(34);
conservation_matrix(35,:) = [1, 0, 0, 1, 1]; %Fe2_Cl_ions = reactants(35);
conservation_matrix(36,:) = [2, 0, 0, 1, 1]; %Fe3_Cl_ions = reactants(36);
conservation_matrix(37,:) = [1, 0, 0, 2, 1];%Fe3_Cl2_ions = reactants(37);

%
options = odeset('RelTol',tol_used,'AbsTol',tol_used);
[times1,eq_concentrations] = ode15s(@water_radiolysis, timeSpan, initialConcentration,...
    options, 0*generationDueToGValues,k);
initialConcentration2 =  eq_concentrations(end,:);
concentrations_elements = [initialConcentration;initialConcentration2]*conservation_matrix
elements_error = concentrations_elements(end,:)-concentrations_elements(1,:)

relative_error = abs(elements_error./concentrations_elements(1,:))
%


[times,concentrations] = ode15s(@water_radiolysis, timeSpan, initialConcentration2,...
    options, generationDueToGValues,k);

end_concentrations(i,:) = interp1(times,concentrations,10^3);
end

run_time = toc
%
% [Fe(III),Cl-, end concentrations]
data_model_2_variable = [trial_map,end_concentrations];
save data_2variable_101by101 data_model_2_variable
%% make concentration maps
for i = 1:37
concentration_maps(:,:,i) = reshape(end_concentrations(:,i),size(Fe_map));
end


%% map out
clear all
load data_2variable_101by101.mat
for i = 1:37
concentration_maps(:,:,i) = reshape(data_model_2_variable(:,i+2),[101,101]);
end
Fe_map = reshape(data_model_2_variable(:,1),[101,101]);
Cl_map = reshape(data_model_2_variable(:,2),[101,101]);
names = ["e^{-}_{aq} (μM)";
    "H_ions";
    "OH_ions";
    "H2O2";
    "HO2_ions";
    "H_mono";
    "OH";
    "O_ions";
    "HO2";
    "O2_ions";
    "O2";
    "H2";
    "O3_ions"
    "O3";
    "HO3";
    "H2O";
    "Cl_ions";
    "ClOH_ions";
    "Cl_mono";
    "Cl2_ions";
    "Cl3_ions";
    "Cl2";
    "Fe2_ions";
    "Fe3_ions";
    "Fe4_ions";
    "Fe5";
    "Fe6_ions";
    "Fe3_OH_ions";
    "Fe3_OH2_ions";
    "Fe3_dimer_ions";
    "Fe2_OH_ions";
    "Fe3_HO2_ions";
    "Fe3_H2O3_ions";
    "O2_2m_ions";
    "Fe2_Cl_ions";
    "Fe3_Cl_ions";
    "Fe3_Cl2_ions";
    ];


%% map all
clf
subplot(6,7,1)
for j = 1:37
    nexttile
plot_species = j;

pcolor(Fe_map,Cl_map,concentration_maps(:,:,plot_species)/10^6)
colormap("jet")
shading interp
box on
ax = gca;
ax.Layer = 'top';
title(names(j),'Interpreter', 'none')
%set(gca,'ColorScale','log')
clim([0,10^round(log10(max(concentration_maps(:,:,plot_species)/10^6,[],'all')),2)])
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
colorbar
xlabel ('Initial Fe(III) (M)')
ylabel ('Initial Cl^{-} (M)')
end



%% Prepping figure titles and cbar scales

names_paper = ["e^{-}_{aq} (μM)";
                "H_2O_2 (mM)";
                "Fe^{2+} (mM)";
                "Cl^{-} (M)";
                "H_2 (M)";
                "O_2 (mM)";
                "Fe^{3+} (mM)";
                "ClOH^- (μM)";
                "H^+ (M)";
                'H^• (μM)';
                'Fe^{4+} (mM)';
                'Cl^• (μM)';
                'HO^- (nM)';
                'HO^• (μM)';
                'Fe^{5+} (mM)';
                'Cl^-_3 (mM)';
                'HO^-_2 (nM)';
                'HO^•_2 (mM)';
                'Fe^{6+} (mM)';
                'Cl_2 (mM)'];
                
conc_exp_paper = [10^-6;
    10^-3;
    10^-3;
    10^0;
    10^0;
    10^-3;
    10^-3;
    10^-6;
    10^0;
    10^-6;
    10^-3;
    10^-6;
    10^-9;
    10^-6;
    10^-3;
    10^-3;
    10^-9;
    10^-3;
    10^-3;
    10^-3];

%% maps for paper

fe_2_map = concentration_maps(:,:,23)+concentration_maps(:,:,35)+concentration_maps(:,:,31);
fe_3_map = concentration_maps(:,:,24)+concentration_maps(:,:,28)+concentration_maps(:,:,29)+...
    concentration_maps(:,:,30)+concentration_maps(:,:,32)+concentration_maps(:,:,33)+...
    concentration_maps(:,:,36)+concentration_maps(:,:,37);

lw = 2;
fz = 10;
figs = 0;
monitor_info = get(0, 'MonitorPositions');
size = 200*[7.5,6];
clf
hold on
t = tiledlayout(4,5,'TileSpacing','compact','Padding','compact');

clim_set = [1*10^-6;
    0.5;
    0.4;
    1;
    1;
    0.4;
    0.5;
    2*10^-6;
    1;
    0.5*10^-4;
    0.2;
    1*10^-5;
    1.5*10^-7;
    0.5*10^-4;
    0.002;
    0.04;
    8.0*10^-8;
    4.0*10^-3;
    0.002;
    2*10^-3]./conc_exp_paper;
names_format = [];

n = 1;

set(gca, 'visible', 'off')
for j = [1,4]
    nexttile
plot_species = j;

p = pcolor(Fe_map,Cl_map,concentration_maps(:,:,plot_species)/10^6/conc_exp_paper(n));
colormap("jet")
shading interp
box on
ax = gca;
ax.Layer = 'top';
title(names_paper(n))
%set(gca,'ColorScale','log')
%clim([0,10^round(log10(max(concentration_maps(:,:,plot_species)/10^6,[],'all')),figs)])
clim([0,clim_set(n)])
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
c = colorbar;
c.Ticks = [0,clim_set(n)/2,clim_set(n)];
%c.Ticks = [0,10^round(log10(max(concentration_maps(:,:,plot_species)/10^6,[],'all')),figs)./2,10^round(log10(max(concentration_maps(:,:,plot_species)/10^6,[],'all')),figs)];
c.TickLength = 0;
xlabel ('Initial Fe(III) (M)')
ylabel ('Initial Cl^{-} (M)')
set(gca,'fontsize',fz,'LineWidth',lw,'FontName','ArialBold','FontWeight','bold')
c.LineWidth = lw;
p = gca;
p.XAxis.TickLength = [0,0];
p.YAxis.TickLength = [0,0];
p.XAxis.TickLabelRotation = 0;

n = n+1;
end

    nexttile
p = pcolor(Fe_map,Cl_map,fe_2_map/10^6/conc_exp_paper(n));
colormap("jet")
shading interp
box on
ax = gca;
ax.Layer = 'top';
title(names_paper(n))
%set(gca,'ColorScale','log')
clim([0,clim_set(n)])
%clim([0,10^round(log10(max(fe_2_map/10^6,[],'all')),figs)])
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
c = colorbar;
c.TickLength = 0;
%c.Ticks = [0,10^round(log10(max(fe_2_map/10^6,[],'all')),figs)./2,10^round(log10(max(fe_2_map/10^6,[],'all')),figs)];
c.Ticks = [0,clim_set(n)/2,clim_set(n)];

xlabel ('Initial Fe(III) (M)')
ylabel ('Initial Cl^{-} (M)')
xticks([0 0.5 1])
yticks([0 0.5 1])
set(gca,'fontsize',fz,'LineWidth',lw,'FontName','ArialBold','FontWeight','bold')
c.LineWidth = lw;
p = gca;
p.XAxis.TickLength = [0,0];
p.YAxis.TickLength = [0,0];
p.XAxis.TickLabelRotation = 0;
n = n+1;
for j = [17,12,11]
    nexttile
plot_species = j;

p = pcolor(Fe_map,Cl_map,concentration_maps(:,:,plot_species)/10^6/conc_exp_paper(n));
colormap("jet")
shading interp
box on
ax = gca;
ax.Layer = 'top';
title(names_paper(n))
%set(gca,'ColorScale','log')
clim([0,clim_set(n)])
%clim([0,10^round(log10(max(concentration_maps(:,:,plot_species)/10^6,[],'all')),figs)])
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
c = colorbar;
c.Ticks = [0,clim_set(n)/2,clim_set(n)];

%c.Ticks = [0,10^round(log10(max(concentration_maps(:,:,plot_species)/10^6,[],'all')),figs)./2,10^round(log10(max(concentration_maps(:,:,plot_species)/10^6,[],'all')),figs)];
c.LineWidth = lw;
c.TickLength = 0;
xlabel ('Initial Fe(III) (M)')
ylabel ('Initial Cl^{-} (M)')
set(gca,'fontsize',fz,'LineWidth',lw,'FontName','ArialBold','FontWeight','bold')
p = gca;
p.XAxis.TickLength = [0,0];
p.YAxis.TickLength = [0,0];
p.XAxis.TickLabelRotation = 0;
n = n+1;
end

    nexttile
p = pcolor(Fe_map,Cl_map,fe_3_map/10^6/conc_exp_paper(n));
colormap("jet")
shading interp
box on
ax = gca;
ax.Layer = 'top';
title(names_paper(n))
%set(gca,'ColorScale','log')
clim([0,clim_set(n)])
%clim([0,10^round(log10(max(fe_3_map/10^6,[],'all')),figs)])
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
c = colorbar;
c.Ticks = [0,clim_set(n)/2,clim_set(n)];

%c.Ticks = [0,10^round(log10(max(fe_3_map/10^6,[],'all')),figs)./2,10^round(log10(max(fe_3_map/10^6,[],'all')),figs)];
xlabel ('Initial Fe(III) (M)')
ylabel ('Initial Cl^{-} (M)')
set(gca,'fontsize',fz,'LineWidth',lw,'FontName','ArialBold','FontWeight','bold')
c.LineWidth = lw;
c.TickLength = 0;
p = gca;
p.XAxis.TickLength = [0,0];
p.YAxis.TickLength = [0,0];
p.XAxis.TickLabelRotation = 0;
n = n+1;

for j = [18,2,6,25,19,3,7,26,21,5,9,27,22]
    nexttile
plot_species = j;

p = pcolor(Fe_map,Cl_map,concentration_maps(:,:,plot_species)/10^6/conc_exp_paper(n));
colormap("jet")
shading interp
box on
ax = gca;
ax.Layer = 'top';
title(names_paper(n))
%set(gca,'ColorScale','log')
clim([0,clim_set(n)])
%clim([0,10^round(log10(max(concentration_maps(:,:,plot_species)/10^6,[],'all')),figs)])
%set(gca, 'XScale', 'log')
%set(gca, 'YScale', 'log')
c = colorbar;
c.TickLength = 0;
c.Ticks = [0,clim_set(n)/2,clim_set(n)];

%c.Ticks = [0,10^round(log10(max(concentration_maps(:,:,plot_species)/10^6,[],'all')),figs)./2,10^round(log10(max(concentration_maps(:,:,plot_species)/10^6,[],'all')),figs)];
xlabel ('Initial Fe(III) (M)')
ylabel ('Initial Cl^{-} (M)')
set(gca,'fontsize',fz,'LineWidth',lw,'FontName','ArialBold','FontWeight','bold')
c.LineWidth = lw;
p = gca;
p.XAxis.TickLength = [0,0];
p.YAxis.TickLength = [0,0];
n = n+1;
end
set(gcf,'color','w')
%set(gcf,'units','pixel','position',[round((monitor_info(3)-size(2))/1.1),0,size(2),size(1)]);

set(gcf,'units','pixel','position',[0,0,900,900]);
%export_fig fe_cl_colormaps_2023_09_25_1round -svg -m3 -painters -nocrop


%% Single map


fe_2_map = concentration_maps(:,:,23)+concentration_maps(:,:,35)+concentration_maps(:,:,31);
fe_3_map = concentration_maps(:,:,24)+concentration_maps(:,:,28)+concentration_maps(:,:,29)+...
    concentration_maps(:,:,30)+concentration_maps(:,:,32)+concentration_maps(:,:,33)+...
    concentration_maps(:,:,36)+concentration_maps(:,:,37);

clim_set = [1*10^-6;
    0.5;
    0.4;
    1;
    1;
    0.4;
    0.5;
    2*10^-6;
    1;
    0.5*10^-4;
    0.2;
    1*10^-5;
    1.5*10^-7;
    0.5*10^-4;
    0.002;
    0.04;
    8.0*10^-8;
    4.0*10^-3;
    0.002;
    2*10^-3]./conc_exp_paper;

paper_figure_index = [1;
    4;
    0;
    17;
    12;
    11;
    0;
    18;
    2;
    6;
    25;
    19;
    3;
    7;
    26;
    21;
    5;
    9;
    27;
    22];

clf
%%
export_fig Fe_vs_Cl_20_Cl2_v1 -svg -m3 -nocrop -painters
%%

for n = 20
if n == 3

     p = pcolor(Fe_map,Cl_map,fe_2_map/10^6/conc_exp_paper(n));

%c_map = sum(conc_map(:,:,2:4),3)
%then plot with p_color
elseif n == 7

    p = pcolor(Fe_map,Cl_map,fe_3_map/10^6/conc_exp_paper(n));


else
plot_species = paper_figure_index(n);
    p = pcolor(Fe_map,Cl_map,concentration_maps(:,:,plot_species)/10^6/conc_exp_paper(n));

end



colormap("jet")

ax = gca;
ax.Layer = 'top';
ax.OuterPosition = [0, 0, 1,1]
ax.InnerPosition = [.14, .14, .8,.78]

title(names_paper(n))
clim([0,clim_set(n)])
c = colorbar;
c.Ticks = [0,clim_set(n)/2,clim_set(n)];
xlabel ('Initial Fe(III) (M)')
ylabel ('Initial Cl^{-} (M)')
c.LineWidth = 4;
p = gca;
p.XAxis.TickLength = [0,0];
p.YAxis.TickLength = [0,0];
yticks([0:.2:1])
xticks([0:.2:1])
set(gcf,'color','w')

box on
set(gca,'Layer','Top','Linewidth', 4,'XMinorTick','Off','YMinorTick','Off')
set(gca,'FontSize',16)
set(gca,'FontName','ArialBold','FontWeight','Bold')
set(gcf,'Renderer','painters')
set(gca,'XColor','black','YColor','black')
shading interp
box(gca,'on');
set(gcf,'color','w')
set(gcf,'Position',[1500 200 600 600])

end

%%

