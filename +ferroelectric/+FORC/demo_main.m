import ferroelectric.FORC.*
clc; clearvars; close all;

% B1500 example
% Open the file for reading, specify file name
raw_file = './FORC_+3V_4m.csv';

area = 50 * 50 * 1e-8;% cm^2
thickness = 9;% nm
% Define the C object
C1 = FORC_single(raw_file,area,thickness,3);

figure;
C1.plot_raw_VIT;
figure;
C1.plot_FORC_trans;


C2 = FORC_single(raw_file,area,thickness,3,'cutRampPoints',3,'Sample_name','test');
figure;
C2.plot_FORC_trans;
figure;
C2.plot_FORC_raw;











