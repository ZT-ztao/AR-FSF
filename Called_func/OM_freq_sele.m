function [A, b, index] = OM_freq_sele(method, sig_sm, sm_snr, threshold_sm, sig_pt, varargin)

if ~isempty(varargin)
    sig_snr = varargin{1};
    threshold_pt = varargin{2};
end


switch method
    case 'Forward'
        noise_band_start_sm = 19e3;
        if_sm_high_snr = [];
        for i=1:3
            temp_sm_snr = sm_snr(:, 26929*(i-1)+1: 26929*i);
            temp_if_sm_high_snr = temp_sm_snr > threshold_sm(i);
            temp_if_sm_high_snr(noise_band_start_sm:end) = 0;
            if_sm_high_snr = [if_sm_high_snr, temp_if_sm_high_snr];
        end
        if_sm_high_snr = logical(if_sm_high_snr);
        index = find(if_sm_high_snr==1);
        A = sig_sm(:,if_sm_high_snr);
        b = sig_pt(if_sm_high_snr);

    case 'FB'
        if_fb_high_snr = [];
        for i=1:3
            temp_sm_snr  = sm_snr(26929*(i-1)+1: 26929*i);
            temp_pt_snr = sig_snr(26929*(i-1)+1: 26929*i);

            temp_if_sm_high_snr = temp_sm_snr  > threshold_sm(i);
            temp_if_pt_high_snr = temp_pt_snr > threshold_pt(i);
            
            temp_if_fb_high_snr  = temp_if_sm_high_snr.' & temp_if_pt_high_snr;
            if_fb_high_snr = [if_fb_high_snr, temp_if_fb_high_snr];
        end
        if_fb_high_snr = logical(if_fb_high_snr);
        index = find(if_fb_high_snr==1);
        A = sig_sm(:,if_fb_high_snr);
        b = sig_pt(if_fb_high_snr);
end


end