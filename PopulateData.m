ecg_sig = [];
result_array = ["temp"];
Fs = 360;
unhealthy_data_struct =  load('data/unhealthy_data.mat');
unhealthy_data = unhealthy_data_struct.unhealthy_data_real;
healthy_data_struct =  load('data/healthy_data.mat');
healthy_data = healthy_data_struct.healthy;

for i=1:24
     if i <10
        strname = strcat('0',int2str(i));
    else
        strname = int2str(i);
     end

    filename = strcat('data/', '1', strname,'m.mat');
    try
        a = load(filename);
    catch ME
        continue
    end

    for ii = 1:length(unhealthy_data(:,1))
        if int2str(unhealthy_data(ii,1)) == strcat('1', strname)
            strt = unhealthy_data(ii,2) - 7200;
            ending = unhealthy_data(ii,2) + 14400;
            ecg_sig = [ecg_sig ;(a.val(1,strt:ending))/200];
            result_array(length(result_array)+1) = "unhealthy";  
        end
    end

    for ii = 1:length(healthy_data(:,1))
        if int2str(healthy_data(ii,1)) == strcat('1', strname)
            strt = healthy_data(ii,2);
            ending = healthy_data(ii,3);
            ecg_sig = [ecg_sig ;(a.val(1,strt:ending))/200];
            result_array(length(result_array)+1) = "healthy";  
        end
    end
end


for i=1:34
    if i <10
        strname = strcat('0',int2str(i));
    else
        strname = int2str(i);
    end

    filename = strcat('data/', '2', strname,'m.mat');
    try
        a = load(filename);
    catch ME
        continue
    end

    for ii = 1:length(unhealthy_data(:,1))
        if int2str(unhealthy_data(ii,1)) == strcat('1', strname)
            strt = unhealthy_data(ii,2) - 7200;
            ending = unhealthy_data(ii,2) + 14400;
            ecg_sig = [ecg_sig ;(a.val(1,strt:ending))/200];
            result_array(length(result_array)+1) = "unhealthy";  
        end
    end

    for ii = 1:length(healthy_data(:,1))
        if int2str(healthy_data(ii,1)) == strcat('1', strname)
            strt = healthy_data(ii,2);
            ending = healthy_data(ii,3);
            ecg_sig = [ecg_sig ;(a.val(1,strt:ending))/200];
            result_array(length(result_array)+1) = "healthy";  
        end
    end
end

t=(0:length(ecg_sig(1,:))-1)/Fs;
plot(t, ecg_sig(1,:));
actual_result = result_array(2:end);

