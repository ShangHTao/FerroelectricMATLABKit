import ferroelectric.PUND.*
clc; clearvars; close all;

% TFA3000 example
% Open the file for reading, specify file name
raw_file_TFA = './PUND_TFA3000.dat';

% Define the TFA_C object
TFA_C = PUND_single_TFA3000(raw_file_TFA);

% plot PUND PV
figure;
TFA_C.plot_PV('black');

% B1500 example
% Open the file for reading, specify file name
raw_file_B1500 = './PUND_B1500.csv';
area = 50 * 50 * 1e-8;% cm^2
thickness = 9;% nm

% Define the B_C object
B_C = PUND_single_B1500(raw_file_B1500,'B1500_test',area,thickness);

% plot PUND PV
figure;
B_C.plot_PV('red');
B_C.get_Pr_Ec











