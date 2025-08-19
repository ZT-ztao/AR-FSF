function adaptive_threshold = cal_OM_apt_threshold(snr,noise_band_start,multiples,type, varargin)
band_one_channel = 26929;

switch length(snr)/band_one_channel
    case 1; channel = 1;
    case 2; channel = 2;
    case 3; channel = 3;
end


for i=1:channel
    snr_channel{i} = snr(band_one_channel*(i-1)+1 : band_one_channel*i);
    noise_band{i} = [noise_band_start: band_one_channel];
end

switch length(varargin)
    case 0
        relax=0.999;
    case 1
        relax=varargin{1};
end



switch type
    case 'sm'
        adaptive_threshold = zeros(channel,1);
        for i=1:channel
            temp_snr = snr_channel{i}(noise_band{i});
            max_noise_value = max(temp_snr);
            relaxed_threshold = relax * max_noise_value;
            adaptive_threshold(i) = multiples * relaxed_threshold;
        end


    case 'pt'
        adaptive_threshold = zeros(channel,1);
        for i=1:channel
            temp_snr = snr_channel{i}(noise_band{i}, :);
            max_noise_value = max(temp_snr);
            relaxed_threshold = relax * max_noise_value;
            adaptive_threshold(i) = multiples * relaxed_threshold;
        end
end




end