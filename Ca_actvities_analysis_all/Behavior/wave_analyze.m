close all; clear all; clc;
FsLabVIEW = 1000;
fish1 = importfile('E:\A_Data_lightsheet\Data_vmat\20190223\fish3\20190223_fish3.txt',1,inf);%,1,1440000
fish1.Properties.VariableNames{3} = 'LabVIEWMeasurement1';
fish1.Properties.VariableNames{5} = 'LabVIEWMeasurement2';
fish1.Properties.VariableNames{7} = 'VarName3';
fish1.Properties.VariableNames{9} = 'LabVIEWMeasurement4';
fish1.Properties.VariableNames{4} = 'CCDSync';
fish1.Properties.VariableNames{6} = 'piezoOut';
fish1.Properties.VariableNames{8} = 'RedLight';
fish1.Properties.VariableNames{10} = 'sCMOSTrigger';
% fish1.Properties.VariableNames{14} = 'US';
figure;plot(fish1.LabVIEWMeasurement*FsLabVIEW,fish1.CCDSync);
figure;plot(fish1.LabVIEWMeasurement,fish1.sCMOSSync);
figure;plot(fish1.LabVIEWMeasurement,fish1.sCMOSTrigger);
% figure;plot(fish1.LabVIEWMeasurement,fish1.US);

