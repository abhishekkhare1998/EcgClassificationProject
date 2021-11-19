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
                ecg = [ecg ;(a.val(1,strt:ending))/200];
                result(length(result)+1) = "unhealthy";  
            end
        end

        for ii = 1:length(healthy(:,1))
            if int2str(healthy(ii,1)) == strcat(series, strname)
                strt = healthy(ii,2);
                ending = healthy(ii,3);
                ecg = [ecg ;(a.val(1,strt:ending))/200];
                result(length(result)+1) = "healthy";  
            end
        end

    end
end