function [snr_sm_mtx, snr_sm_vct] = cal_OM_snr_SM(S_sig, S_bg, type, varargin)
dim_sig = size(S_sig);
dim_bg  = size(S_bg);
S_sig = reshape(S_sig, dim_sig(1), dim_sig(2)*dim_sig(3));
S_bg = reshape(S_bg, dim_bg(1), dim_bg(2)*dim_bg(3));


%% type of snr calculations
[OpenMPI_v_ori, ~] = cal_OpenMPI_velocity;
OpenMPI_v_norm = OpenMPI_v_ori/max(OpenMPI_v_ori(:));

if isempty(varargin)
    tau_v = 0.9;
else 
    tau_v = varargin{1};
end


if_high_v_pos  = OpenMPI_v_norm(:) > tau_v*max(OpenMPI_v_norm(:));
high_v_pos = find( if_high_v_pos == 1);

S_sig_high_v = S_sig(high_v_pos, :);
S_bg_high_v = S_bg(high_v_pos, :);


switch type
    case 'pale'
        fprintf('---------- calcalating snr_sm pale ---------- \n')
        snr_sm_mtx = [];
        snr_sm_vct = sqrt(sum(abs(S_sig).^2,1))./sqrt(sum(abs(S_bg).^2,1));
    case 'vc-pale'
        fprintf('---------- calcalating snr_sm vc-pale ---------- \n')
        snr_sm_mtx = abs(S_sig_high_v)./sqrt(sum(abs(S_bg_high_v).^2,1)/size(S_bg_high_v,1));
        snr_sm_vct = sum((OpenMPI_v_norm(high_v_pos) .* snr_sm_mtx).^2, 1) / sum(if_high_v_pos);
    case 'std'
        fprintf('---------- calcalating snr_sm std ---------- \n')
        snr_sm_mtx = [];
        snr_sm_vct = sqrt(sum(abs(S_sig-S_bg).^2,1))./(sqrt(dim_bg(1))*std(S_bg));
    case 'vc-std'
        fprintf('---------- calcalating snr_sm vc-std ---------- \n')
        snr_sm_mtx = abs(S_sig_high_v - S_bg_high_v) ./ std(S_bg_high_v);
        snr_sm_vct = sum((OpenMPI_v_norm(high_v_pos) .* snr_sm_mtx).^2, 1) / sum(if_high_v_pos);
end


end

