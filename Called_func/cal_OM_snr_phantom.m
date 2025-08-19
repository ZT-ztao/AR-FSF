function [snr_pt_vct] = cal_OM_snr_phantom(U_freq_pt, U_freq_bg, type)

dim_sig = size(U_freq_pt);
dim_bg  = size(U_freq_bg);

pt_sig = reshape(U_freq_pt, dim_sig(1)*dim_sig(2), dim_sig(3));
pt_bg = reshape(U_freq_bg, dim_bg(1)*dim_bg(2), dim_bg(3));

switch type
    case 'pale'
        snr_pt_vct = sqrt(sum(abs(pt_sig).^2,2))./sqrt(sum(abs(pt_bg).^2,2));
    case 'std'
        snr_pt_vct = sqrt(sum(abs(pt_sig-pt_bg).^2,2))./(std(pt_bg,0,2)*sqrt(dim_bg(3)));
    case 'vc-pale'
        snr_pt_vct = (sqrt(sum(abs(pt_sig).^2,2))./sqrt(sum(abs(pt_bg).^2,2))).^2;
    case 'vc-std'
        snr_pt_vct = (sqrt(sum(abs(pt_sig-pt_bg).^2,2))./(std(pt_bg,0,2)*sqrt(dim_bg(3)))).^2;
end


end
