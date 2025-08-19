%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an example implementation of the Adaptive and Robust Frequency Selection Framewrok (AR-FSF)
% In this example, used data are from Open-MPI dataset
% ----Example data can download from https://media.tuhh.de/ibi/openMPIData/data/calibrations/
% ----Detailed data explanations can be found in DOI:10.1016/j.dib.2019.104971
% In this example, data need to be previously downloaded to local files, and then define the (file_path / used_data_num) variables 

clear
%% imprt data
% loading the system matrix data and the phantom data in the Open-MPI dataset. 
fprintf('.....................Loading data.....................\n')

file_path_sm = '..........';  % the file path of the system matrix data
file_path_pt = '..........';  % the file path of the phantom data
used_data_num = 3;            % the used data number

[S_meas, S_bg, ~, ~, ~] = Read_OpenMPI_SM(file_path_sm, used_data_num);
[U_reso_sig, U_reso_bg, ~, ~]= Read_OpenMPI_phantom('resolution', used_data_num);
fprintf('.....................Loading data Finish!.....................\n')

%% calculating SM-SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the function cal_OM_snr_SM, there are four options fpr calculating the snr for system matrix. 
% ----------pale   : SNR_{nrg} feature;
% ----------std    : SNR_{std} feature;
% ----------vc-pale: velocity-corrected SNR_{nrg} feature;
% ----------vc-std : velocity-corrected SNR_{std} feature;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('.....................calculating SM-SNR........................\n')
[~, snr_sm_vcpale_vct]  = cal_OM_snr_SM(S_meas,S_bg,'vc-pale');

%% calculating Pt-SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the function cal_OM_snr_SM, there are four options fpr calculating the snr for phantom. 
% ----------pale   : SNR_{nrg} feature;
% ----------std    : SNR_{std} feature;
% ----------vc-pale: velocity-corrected SNR_{nrg} feature;
% ----------vc-std : velocity-corrected SNR_{std} feature;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('.....................calculating Pt-SNR........................\n')
snr_pt_vcpale_vct  = cal_OM_snr_phantom(U_reso_sig, U_reso_bg,'vc-pale');


%% calculating selection threshold of the feature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  theta can be defined according to the imaging requirement :
%  high-resolution imaging        --> a high value above the noise level;
%  high-speed / low-noise imaging --> a low value above the noise level;
%  if there are no specific requirement, theta = 1 are be used. 
%  start_fre : the start frequency of the used noise band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = 1; 
% theta = 2; 
% theta = 3; 

relax = 0.9;
% relax = 0.99;
% relax = 0.999;

start_fre = [19e3, 19e3, 19e3];

apt_thre_sm_vcpale   = cal_OM_apt_threshold(snr_sm_vcpale_vct, start_fre, theta, 'sm', relax);
apt_thre_pt_vcpale   = cal_OM_apt_threshold(snr_pt_vcpale_vct, start_fre, theta, 'pt', relax);

%% selecting frequency components for reconstruction
SM = S_meas-S_bg;
U = mean(U_reso_sig-U_reso_bg,3);
[A_vcpale,b_vcpale,idx_vcpale] = OM_freq_sele('FB', SM, snr_sm_vcpale_vct, apt_thre_sm_vcpale, U, snr_pt_vcpale_vct, apt_thre_pt_vcpale);

%% Image reconstruction
% use the Kaczmarz method DOI: 10.1088/0031-9155/55/6/003
lambd = 1e-3;
iter = 20;
[x_pale, time_pale]    = OpenMPI_kz(A_vcpale,  b_vcpale,  lambd, iter);

