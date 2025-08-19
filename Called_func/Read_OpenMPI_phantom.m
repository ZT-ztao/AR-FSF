function [U_freq_pt, U_freq_bg, U_freq_pt_mean, U_freq_bg_mean] = Read_OpenMPI_phantom(file_path, phantom, read_num)

file_name = [file_path, phantom, '\', num2str(read_num),'.mdf'];
isBG_U = h5read(file_name, '/measurement/isBackgroundFrame'); 
U_time = h5read(file_name, '/measurement/data');    
U_time_pt = squeeze(U_time(:,:,:,isBG_U == 0));
U_time_bg = squeeze(U_time(:,:,:,isBG_U == 1));
U_time_pt_mean = mean(U_time_pt,3);
U_time_bg_mean = mean(U_time_bg,3);

U_freq_pt = fft(cast(U_time_pt,'double')); 
U_freq_pt = U_freq_pt(1:(size(U_freq_pt,1)/2+1),:,:); 
U_freq_bg = fft(cast(U_time_bg,'double'));  
U_freq_bg = U_freq_bg(1:(size(U_freq_bg,1)/2+1),:,:); 


U_freq_pt_mean = fft(cast(U_time_pt_mean,'double'));
U_freq_pt_mean = U_freq_pt_mean(1:(size(U_freq_pt_mean,1)/2+1),:);
U_freq_bg_mean = fft(cast(U_time_bg_mean,'double'));  
U_freq_bg_mean = U_freq_bg_mean(1:(size(U_freq_bg_mean,1)/2+1),:);


end