function [S_meas, S_bg, snr, isBG, dim] = Read_OpenMPI_SM(file_path, read_num)

file_name = [file_path, num2str(read_num),'.mdf'];
S_ri = h5read(file_name, '/measurement/data');       
S_Ori = double(complex(S_ri.r,S_ri.i));              

isBG = h5read(file_name, '/measurement/isBackgroundFrame');
dim = h5read(file_name, '/calibration/size');               
snr = h5read(file_name, '/calibration/snr');             


S_meas= S_Ori(isBG == 0,:,:); 
S_bg = S_Ori(isBG == 1,:,:); 
end