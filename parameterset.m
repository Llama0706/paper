
time = acquired.data.time;

stimulus_info_1 = acquired.stim(1,1).data;%Exporting stimulus information
stimulus_info_2 = acquired.stim(1,2).data;
stimulus_info_3 = acquired.stim(1,3).data;
stimulus_sequence_1 = zeros(length(time), 1);
stimulus_sequence_2 = zeros(length(time), 1);
stimulus_sequence_3 = zeros(length(time), 1);
for i = 1:size(stimulus_info_1, 1)
    start_time = stimulus_info_1(i, 1); 
    duration = stimulus_info_1(i, 2);  
    end_time = start_time + duration; 
    
    % Find the index corresponding to the start time and end time
    start_idx = find(time == start_time, 1);
    end_idx = find(time == end_time, 1);
    
    if isempty(start_idx) || isempty(end_idx)
        disp(['Warning: Start or end time not found in time series mapping']);
        continue;
    end
    
    % Setting the stimulus sequence to 1 in the index range between start time and end time
    stimulus_sequence_1(start_idx:end_idx) = 1;
end%Populate the content of the newly created stimulus sequence according to the stimulus sequence information, setting the value of the stimulus position to 1

%Generating a second stimulus sequence
for i = 1:size(stimulus_info_2, 1)
    start_time = stimulus_info_2(i, 1); 
    duration = stimulus_info_2(i, 2);  
    end_time = start_time + duration;
    
   
    start_idx = find(time == start_time, 1);
    end_idx = find(time == end_time, 1);
    
    if isempty(start_idx) || isempty(end_idx)
        disp(['Warning: Start or end time not found in time series mapping']);
        continue;
    end
    
    stimulus_sequence_2(start_idx:end_idx) = 1;
end



for i = 1:size(stimulus_info_3, 1)
    start_time = stimulus_info_3(i, 1);
    duration = stimulus_info_3(i, 2);   
    end_time = start_time + duration; 
    
  
    start_idx = find(time == start_time, 1);
    end_idx = find(time == end_time, 1);
    
    if isempty(start_idx) || isempty(end_idx)
        disp(['Warning: Start or end time not found in time series mapping']);
        continue;
    end
    
    stimulus_sequence_3(start_idx:end_idx) = 1;
end

hbO_data = outResults.dc.dataTimeSeries(:, 1:61);
max_duration = min(30, max(time));
% Defining parameter ranges
tau_p_values = 1:20;
tau_d_values = 3:10;

% Initialize the maximum average t-value and the corresponding parameter
max_avg_t = -inf;
best_tau_p = 0;
best_tau_d = 0;
t = time(time <= max_duration);
A = 4;
% Loop over all combinations of tau_p and tau_d
for tau_p = tau_p_values
    for tau_d = tau_d_values
        % Recalculation of HRF
        hrf = (t.^tau_p).*exp(-t)/factorial(tau_p) - (t.^(tau_p + tau_d).*exp(-t))/(A*factorial(tau_p + tau_d));
        
       % Calculate the convolution result
convolution_result_1 = conv(hrf, stimulus_sequence_1);
convolution_result_2 = conv(hrf, stimulus_sequence_2);
convolution_result_3 = conv(hrf, stimulus_sequence_3);
num_channels = 61;
betas_hbO_S1 = zeros(num_channels, 2);
betas_hbO_S2 = zeros(num_channels, 2);
betas_hbO_S3 = zeros(num_channels, 2);

for i = 1:num_channels
    X_HbO2_1 = [convolution_result_1, ones(length(convolution_result_1), 1)];
    X_HbO2_2 = [convolution_result_2, ones(length(convolution_result_2), 1)];
    X_HbO2_3 = [convolution_result_3, ones(length(convolution_result_3), 1)];
    X_HbO2_1 = X_HbO2_1(1:length(hbO_data), :);
    X_HbO2_2 = X_HbO2_2(1:length(hbO_data), :);
    X_HbO2_3 = X_HbO2_3(1:length(hbO_data), :);
    betas_hbO_S1(i, :) = X_HbO2_1 \ hbO_data(:, i); 
    betas_hbO_S2(i, :) = X_HbO2_2 \ hbO_data(:, i);
    betas_hbO_S3(i, :) = X_HbO2_3 \ hbO_data(:, i);
end

        
        beta_values1 = betas_hbO_S1(:, 1);
        beta_values2 = betas_hbO_S2(:, 1);
        beta_values3 = betas_hbO_S3(:, 1);
        
zero_indices = find(plotOptions.blockActive{1, 1}(:, 3) == 0);

for i = 1:numel(zero_indices)
    if zero_indices(i) > 61
        zero_indices(i) = zero_indices(i) - 61;
    end
end


beta_values1(zero_indices) = [];
beta_values2(zero_indices) = [];
beta_values3(zero_indices) = [];


[h1, p1, ~, stats1] = ttest(beta_values1, 0, 'Tail', 'right', 'Alpha', 0.05);
[h2, p2, ~, stats2] = ttest(beta_values2, 0, 'Tail', 'right', 'Alpha', 0.05);
[h3, p3, ~, stats3] = ttest(beta_values3, 0, 'Tail', 'right', 'Alpha', 0.05);

        
    
        avg_t = (stats1.tstat + stats2.tstat + stats3.tstat) / 3;
        
        if avg_t > max_avg_t
            max_avg_t = avg_t;
            best_tau_p = tau_p;
            best_tau_d = tau_d;
        end
    end
end

disp(['最佳参数 tau_p: ', num2str(best_tau_p)]);
disp(['最佳参数 tau_d: ', num2str(best_tau_d)]);
disp(['最大平均 t 值: ', num2str(max_avg_t)]);