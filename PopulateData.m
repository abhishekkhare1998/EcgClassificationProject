ecg_sig = [];
result_array = ["temp"];
Fs = 360;
unhealthy_data_struct =  load('data/unhealthy_data.mat');
unhealthy_data = unhealthy_data_struct.unhealthy_data_real;
healthy_data_struct =  load('data/healthy_data.mat');
healthy_data = healthy_data_struct.healthy;


[ecg_sig, result_array] = populate(ecg_sig,result_array,'1',24,unhealthy_data, healthy_data);
[ecg_sig, result_array] = populate(ecg_sig,result_array,'2',34,unhealthy_data, healthy_data);


t=(0:length(ecg_sig(1,:))-1)/Fs;
plot(t, ecg_sig(1,:));
actual_result = result_array(2:end);

energy_mat = [];
for i=1:length(actual_result)
    [wavelet_energy] = compute_dwt('db2', 10, ecg_sig(i,:));
    energy_mat = [energy_mat;transpose(wavelet_energy)];
end

%plot_all(ecg_sig)

function [] = plot_all(signal_matrix)
    for i = 1:73
        figure;
        plot(signal_matrix(i,:))
    end
end

function [processed_signal] = pre_processing(signal)
%[processed_signal] = remove_mean(signal);
%[processed_signal] = remove_baseline(signal);
%[processed_signal] = lowpassfilter(removed_mean);
processed_signal = signal;
end

function [filtered_signal] = lowpassfilter(signal)
lpf = load('filters/LP20.mat');
filtered_signal = conv(signal, lpf.Hlp.Numerator);
end

function [removed_avg] = remove_mean(signal)
mean = sum(signal)/length(signal);
for i = 1:length(signal)
    removed_avg(i) = signal(i) - mean;
end
end

function [removed_base] = remove_baseline(signal)
removed_base = [];
i = 1;
incr = 300;
    while(i < length(signal))
        if i + incr > length(signal)
            incr = length(signal) - i;
        end
        local_mean = sum(signal(i:i+incr))/(incr+1);
        for j = i:i+incr
            removed_base(j) = signal(j) - local_mean;
        end
        i = i + incr;
    end
end

function [ecg, result] = populate(ecg, result, series, iter, ...
    unhealthy, healthy)
    for i=1:iter
        if i < 10
            strname = strcat('0',int2str(i));
        else
            strname = int2str(i);
        end

        filename = strcat('data/', series, strname,'m.mat');
        try
            a = load(filename);
        catch ME
            continue
        end

        for ii = 1:length(unhealthy(:,1))
            if int2str(unhealthy(ii,1)) == strcat(series, strname)
                strt = unhealthy(ii,2) - 7200;
                ending = unhealthy(ii,2) + 14400;
                %[signal_of_interest] = remove_mean(a.val(1,strt:ending));
                [signal_of_interest] = pre_processing(a.val(1,strt:ending));
                ecg = [ecg ; (signal_of_interest)/200];

                result(length(result)+1) = "unhealthy";  
            end
        end

        for ii = 1:length(healthy(:,1))
            if int2str(healthy(ii,1)) == strcat(series, strname)
                strt = healthy(ii,2);
                ending = healthy(ii,3);
                [signal_of_interest] = pre_processing(a.val(1,strt:ending));
                ecg = [ecg ;(signal_of_interest)/200];
                result(length(result)+1) = "healthy";  
            end
        end

    end
end

function [wavelet_energy] = compute_dwt(type, level, signal)
    wavelet_energy = [];
    [coefficients,len] = wavedec(signal,level,type);
    j = 1;
    for i = 1:length(len)-1
        wavelet_energy = [wavelet_energy; sum(coefficients(j:j+len(i)-1).^2)/(len(i))];
    end

end
