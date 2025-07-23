classdef photocurrent_single
    properties
        CVS_file_path = '';
        Setup_title = '';
        Voltage_raw = [];
        Time_raw = [];
        Current_raw = [];

        Current_mean = 0;
        Current_absmax = 0;
    end

    methods
        function obj = photocurrent_single(cvs_file_path)
            obj.CVS_file_path = cvs_file_path;
            obj = extract_raw_data(obj);
            obj = cal_mean_and_absmax_current(obj);

        end

        function obj = extract_raw_data(obj)
            file_path = obj.CVS_file_path;

            fileID = fopen(file_path, 'r');
            V = [];
            I = [];
            T = [];
            line = fgets(fileID);
            while ischar(line)
                % N_line = N_line + 1;
                if startsWith(line, 'SetupTitle, WGFMU Pattern Editor Child DataDisplay')
                    break;  % 结束读取后面重复数据
                elseif startsWith(line, 'SetupTitle')
                    splitLine = split(line, ',');
                    if numel(splitLine) >= 2
                        setup_title = strtrim(splitLine{2});
                    end
                elseif startsWith(line, 'DataValue')
                    splitLine = split(line, ',');
                    V(end + 1) = str2double(splitLine{2});
                    I(end + 1) = str2double(splitLine{4});
                    T(end + 1) = str2double(splitLine{3});
                    % raw_data{end+1} = splitLine;  % 存储为 cell 的 cell，每行一个单元
                end
                line = fgetl(fileID);                % 读取下一行
            end
            fclose(fileID);
            obj.Setup_title = setup_title;
            obj.Voltage_raw = V.';
            obj.Current_raw = I.';
            obj.Time_raw = T.';
        end

        function obj = cal_mean_and_absmax_current(obj)
            I = obj.Current_raw;

            obj.Current_mean = mean(I);
            obj.Current_absmax = max(abs(I));
        end
        function obj =  plot_raw_IT(obj,filter_window)

            set_title = obj.Setup_title;
            t_s = obj.Time_raw;
            % v = obj.Voltage_raw;
            c_pA = obj.Current_raw*1e12;
            if ~exist('filter_window','var')
                filter_window = round(0.1*numel(c_pA));
            end


            plot(t_s, c_pA,'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'none', 'MarkerSize', 8, 'Color', 'blue');
            hold on;

            I_filtered = movmean(c_pA, filter_window);

            plot(t_s, I_filtered,'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'none', 'MarkerSize', 10, 'Color', 'red');
            title(set_title,FontSize=12);
            xlabel('Time (s)');
            ylabel('Current (pA)');
            set(gca, 'FontSize', 12);  % Set font size
            legend('raw data', 'filter data', 'FontSize',12,'Location','Best');
        end
    end
end