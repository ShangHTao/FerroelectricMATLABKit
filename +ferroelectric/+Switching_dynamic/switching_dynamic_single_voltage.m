classdef switching_dynamic_single_voltage
    % 测试单个测试电压的switching dynamic数据处理类
    properties
        N_test = 0;

        Time = {};% unit s 含有1xN_test个原胞，每个原胞包含了时间的向量
        Voltage = {};% unit V
        Current = {};% unit A
        Setup_Title = {};% 提取raw data种的setup title

        Time_split = {};% unit s
        Voltage_split = {};% unit V
        Current_split = {};% unit A

        Current_subtracted = {};% unit A 第一段电流减第二段电流
        Polarity = 0;% -1为负电压，1为正电压
        Pulse_time = [];
        P_remnant = [];


        Fitted_t0 = 0;%拟合的数据
        Fitted_w = 0;%拟合的数据
        Fitted_pt = [];%拟合的Pulse_time
        Fitted_P = [];

        Thickness = 0;% unit nm
        % 通常将测试电压作为文件夹归类第一个电容文件夹C1，测试电压2.5V  .../C1/2.5V
        Label_device = '';% e.g. 'C1'
        Label_voltage  = '';% e.g. '2.5V'
        Label_pulse_time = {};% e.g. {'10ns'} {}

        Folder_path = '';
        Area = 0;
        V_bias = 0;
    end

    properties (Access = private)
        Voltage_th = 0.1;
        Window = 5;
    end


    methods
        function obj = switching_dynamic_single_voltage(Folder_path,area,thickness,v_bias)
            % folderPath e.g.'C:\Users\shanght\Desktop\Fe\experimental_data\radiation_SD\He1e15\C1\1V'
            % area unit cm2
            % V_bias unit V
            % thickness unit nm thickness = 10nm
            obj.Folder_path = Folder_path;
            obj.Area = area;
            obj.V_bias = v_bias;
            obj.Thickness = thickness;
            % 使用 filesep 分隔路径
            parts = strsplit(Folder_path, filesep);
            obj.Label_device = parts{end-1};   % 倒数第二项 -> e.g. C1
            obj.Label_voltage = parts{end};    % 倒数第一项  -> e.g. 1V

            obj = extract_raw_data(obj);
            obj = split_data(obj);
            obj = subtracted_current(obj);
            obj = Q_integration(obj);
            obj = auto_NLS_fitting(obj,2);
        end

        function obj = extract_raw_data(obj)
            data_path = obj.Folder_path;
            fileList = dir(fullfile(data_path, 'Read*.csv')); % all files
            N_fileList = numel(fileList);
            condition_pulse_time = cell(N_fileList,1);
            pulse_t = zeros(N_fileList,1);
            % 按照脉冲时间重新排序
            for i = 1:N_fileList
                filename = fileList(i).name;
                str_match = regexp(filename, '(?<=Force ).*?(?= )', 'match');
                if isempty(str_match)
                    error('未找到匹配字段: %s', filename);
                end
                str = str_match{1};
                condition_pulse_time{i} = str;
                if contains(str, 'ns')
                    value = str2double(extractBefore(str, 'ns')) * 1e-9;
                elseif contains(str, 'us')
                    value = str2double(extractBefore(str, 'us')) * 1e-6;
                elseif contains(str, 'ms')
                    value = str2double(extractBefore(str, 'ms')) * 1e-3;
                else
                    value = NaN; % Handle unsupported units or invalid format
                end
                pulse_t(i) = value;% Unsupported

            end
            [~, sortIdx] = sort(pulse_t);% 按照脉冲值排序
            pulse_t = pulse_t(sortIdx);
            fileList = fileList(sortIdx);
            condition_pulse_time = condition_pulse_time(sortIdx);

            time = cell(N_fileList,1);
            voltage = cell(N_fileList,1);
            current = cell(N_fileList,1);
            setup_title = cell(N_fileList,1);
            for i = 1:N_fileList
                filename = fileList(i).name;% Get the current file name
                full_path = [data_path filesep filename];
                fileID = fopen(full_path, 'r');
                row_number = 0; %row number ot each *.dat file
                Row_sum = {};
                title = {}; % extract recipe name, such as "I-V 5us 100kHz 3V - advanced"
                TIV_data = []; % All time, Voltage, Current data
                line = fgetl(fileID); % Read each line using fgetl until the end of the file
                while ischar(line)
                    if strncmp(line, 'DataValue', numel('DataValue'))
                        splitLine2 = split(line, ',');
                        TIV_data = [TIV_data, splitLine2];
                        row_number = row_number + 1;
                    elseif strncmp(line, 'SetupTitle', numel('SetupTitle'))
                        % Extract the next text after "SetupTitle"
                        splitLine1 = split(line, ',');
                        newtext = strtrim(splitLine1{2});
                        title = [title, newtext];
                        Row_sum = [Row_sum, row_number];
                        row_number = 0;
                    end
                    % Read the next line
                    line = fgetl(fileID);
                end
                fclose(fileID);
                % keep only the even-indexed elements (2nd, 4th, 6th, etc.)
                Row_sum = cell2mat(Row_sum(2:2:end));
                % Extend the columns by repeating each element
                Row_sum = reshape(repelem(Row_sum, 2), 1, []);
                %Transpose the array
                TIV_data = TIV_data';
                % Delete the 1st,5th,6th columns
                TIV_data(:, [1, 5]) = [];
                TIV_data = TIV_data(1:row_number,:);
                % Move the original 4th column to the 1st position
                TIV_data = [TIV_data(:, 3), TIV_data(:, 1:2)];
                TIV_data = cellfun(@str2double, TIV_data);% convert sring to double
                TIV_data(:, 3) = TIV_data(:, 3)*-1;% WGFMU2 unit detected *-1 current !!!!

                time{i} = TIV_data(:,1);%1x3 cell, each element is a 3x1000 cell
                voltage{i} = TIV_data(:,2);
                current{i} = TIV_data(:,3);
                setup_title{i} = title{1};
            end
            obj.N_test = N_fileList;
            obj.Time = time;
            obj.Voltage = voltage;
            obj.Current = current;
            obj.Setup_Title = setup_title;
            obj.Pulse_time = pulse_t;
            obj.Label_pulse_time = condition_pulse_time;
        end

        function obj = split_data(obj,voltage_th,window)
            if nargin >= 2
                obj.Voltage_th = voltage_th;
            end
            if nargin >= 3
                obj.Window = window;
            end
            voltage_th = obj.Voltage_th;
            window = obj.Window;
            N = obj.N_test;
            v_bias = obj.V_bias;
            raw_data = [obj.Time, obj.Voltage, obj.Current];

            time_split = cell(N,1);
            voltage_split = cell(N,1);
            current_split = cell(N,1);
            polarity = zeros(N,1);
            for i = 1:N
                [time,voltage,current] = raw_data{i,:};


                if abs(min(voltage)) > 0.5
                    polarity(i) = -1;
                elseif abs(max(voltage)) > 0.5
                    polarity(i) = 1;
                end
                voltage = voltage+v_bias;
                v = voltage;
                dv = abs(diff(v));
                ddv = diff(dv);

                [~, idx] = max(ddv(window:end-window));
                split_idx = idx + window;
                if abs(voltage(split_idx))<abs(voltage_th+v_bias)
                    % disp('Triangular wave separation successful');
                else
                    error('Triangular wave separation failed! Adjast voltage threshold');
                end
                %鲁棒性未优化，必须保证前一段电压持续时间大于后一段
                time_split{i}(:,1) = time(1:split_idx-1,:);%时间划分第一列
                time_split{i}(:,2) = time(split_idx:2*(split_idx-1),:);%时间划分第二列
                voltage_split{i}(:,1) = voltage(1:split_idx-1,:);
                voltage_split{i}(:,2) = voltage(split_idx:2*(split_idx-1),:);
                current_split{i}(:,1) = current(1:split_idx-1,:);
                current_split{i}(:,2) = current(split_idx:2*(split_idx-1),:);
            end
            obj.Time_split = time_split;
            obj.Voltage_split = voltage_split;
            obj.Current_split = current_split;
            obj.Polarity = polarity;
        end

        function obj = subtracted_current(obj)
            N = obj.N_test;
            current_split = obj.Current_split;
            substructed_current = cell(N,1);
            for i = 1:N
                % current P U N D
                current_PN = current_split{i}(:,1);
                current_UD = current_split{i}(:,2);
                substructed_current{i} = current_PN - current_UD;
            end
            obj.Current_subtracted = substructed_current;
        end

        function plot_subtracted_IV(obj,sequence)
            % figure;
            times = obj.Time_split{sequence}(:,1);
            voltage = obj.Voltage_split{sequence}(:,1);
            current_1 = obj.Current_split{sequence}(:,1);
            current_2 = obj.Current_split{sequence}(:,2);
            current_3 = obj.Current_subtracted{sequence};
            fullTitle = [obj.Label_device,'_',obj.Label_voltage,'_',obj.Label_pulse_time{sequence}];

            plot(times, current_1,'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'none', 'MarkerSize', 8, 'Color', 'Blue');
            hold on;
            plot(times, current_2,'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'none', 'MarkerSize', 8, 'Color', 'Green');
            hold on;
            plot(times, current_3,'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'none', 'MarkerSize', 8, 'Color', [0.6, 0.4, 0.2]);
            yline(0,'LineWidth', 1.5,'Color', 'black'); % y = 0
            xline(0,'LineWidth', 1.5,'Color', 'black'); % x = 0
            xlabel('Time (s)', 'FontSize',16);
            ylabel('Current (A)', 'FontSize',16,'Color', 'black');
            % legend('1st Pulse','2nd Pulse','Subtraction','Voltage', 'FontSize',16, 'Location','northwest');
            % Place the title with appropriate font size using the concatenated string
            title(fullTitle, 'FontSize', 16,Interpreter='none'); % Show title at the top of subplots

            ylimit = max(abs(ylim));
            ylim([-ylimit ylimit]);
            xlimit = max(abs(xlim));
            if rem(xlimit,2) == 1
                xlimit = xlimit + 1;
            end
            xlim([0 xlimit]);
            set(gca, 'FontSize', 16);  % Set font size
            grid on;
            grid minor;  % Use minor grid lines
            set(gca, 'LineWidth', 1.5, 'GridColor', 'black');  % Set line thickness and color
            % Plot Voltage signal on the right side axis
            hold on; % Enable holding the plot
            % Plot the second column (data(:,2)) against the first column (data(:,1)) on the right y-axis
            yyaxis right;
            plot(times, voltage,'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'none', 'MarkerSize', 8, 'Color', 'Black');
            if mean(voltage) < 0
                legend('1st Pulse','2nd Pulse','Subtraction','Voltage', 'FontSize',16, 'Location','northwest');
            elseif mean(voltage) > 0
                legend('1st Pulse','2nd Pulse','Subtraction','Voltage', 'FontSize',16, 'Location','southwest');
            end
            % Set appropriate label for the right y-axis plot
            ylabel('Voltage (V)', 'FontSize',16,'Color', 'Black');
            ylim([-4 4]);
        end

        function obj = Q_integration(obj)
            N = obj.N_test;
            area = obj.Area;

            times_split = obj.Time_split;
            current_subtracted = obj.Current_subtracted;

            Qr = zeros(N,1);
            for i = 1:N
                times = times_split{i}(:,1);
                current = current_subtracted{i};
                Qr(i) = trapz(times, current);% unit：库仑（C）
            end
            obj.P_remnant =  Qr*1e6/area;
        end

        function plot_P_remnant(obj)
            % figure;
            pulse_time = obj.Pulse_time;
            p_remnant = obj.P_remnant;
            fullTitle = [obj.Label_device,'_',obj.Label_voltage];

            semilogx(pulse_time, p_remnant, 'o', 'MarkerSize', 8, 'Color', 'black');
            hold on;

            title(fullTitle,'FontSize', 16,'Color', 'black',Interpreter='none');
            yline(0, 'LineWidth', 1.5, 'Color', 'black'); % y = 0
            xline(0, 'LineWidth', 1.5, 'Color', 'black'); % x = 0
            xlabel('Pulse Duration (s)', 'FontSize', 16);
            ylabel('$Q_r$ ($\mu C/cm^2$)', 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'black');
            ylim([-60 60]);
            xlim([1e-8 1e-2]);
            grid on;
            grid minor; % Use minor grid lines
            set(gca, 'LineWidth', 1.5, 'GridColor', 'black'); % Set line thickness and grid color
        end

        function obj = auto_NLS_fitting(obj,n,span)
            if nargin <3
                span=5;% 用于计算经过平滑处理的值的数据点数
            end
            polarity = obj.Polarity(1);
            t = obj.Pulse_time;
            p = obj.P_remnant;
            p_max = max(abs(p)); %fp - > fraction polarity
            p_norm = abs(p/p_max);
            fp_smooth = smooth(p_norm, span, 'moving');
            dpdlogt = gradient(fp_smooth)./gradient(log10(t));
            [~,index] = max(abs(dpdlogt));
            t0_init = t(index);
            %求初始化w0
            % k = dpdlogt(index),logt0 = log10(t0_init);
            % fp = fp_smooth(index)
            % y - fp = k * (logt - logt0)
            % y -> 0 求出与logt轴交点 logt1 = -fp/k + logt0
            % w = logt0 - logt1 = fp/k
            k = dpdlogt(index);
            fp = fp_smooth(index);
            w0_init = abs(fp/k);

            initial_params = [w0_init,t0_init];

            % chi_squared = chi_squared_cal(t,p_norm,n,w0_init,t0_init)
            chi_squared = @(params) sum((p_norm - switching_dynamic_single_voltage.NLS_model(t,n,params(1),params(2))).^2); % Define the chi-squared function

            fitted_params = fminsearch(chi_squared, initial_params);

            t_fit = logspace(log10(t(1)), log10(t(end)), 1000);
            fraction_P_fit = switching_dynamic_single_voltage.NLS_model(t_fit,n, fitted_params(1), fitted_params(2));
            p_fit = polarity * p_max * fraction_P_fit;

            obj.Fitted_w = fitted_params(1);
            obj.Fitted_t0 = fitted_params(2);
            obj.Fitted_pt = t_fit;
            obj.Fitted_P = p_fit;
            % semilogx(t_fit,p_fit);
            % hold on;
            % semilogx(t,p);
        end

        function plot_fitted_P_remnant(obj)
            fitted_pt = obj.Fitted_pt;
            fitted_P = obj.Fitted_P;
            semilogx(fitted_pt, fitted_P, 'LineWidth',1,'LineStyle','-', 'Color','black');
            fultitle = [obj.Label_device '_' obj.Label_voltage ' fitted data'];
            title(fultitle,'FontSize', 16,'Color', 'black',Interpreter='none');
            yline(0, 'LineWidth', 1.5, 'Color', 'black'); % y = 0
            xline(0, 'LineWidth', 1.5, 'Color', 'black'); % x = 0
            xlabel('Pulse Duration (s)', 'FontSize', 16);
            ylabel('$Q_r$ ($\mu C/cm^2$)', 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'black');
            ylim([-60 60]);
            xlim([1e-8 1e-2]);
            grid on;
            grid minor; % Use minor grid lines
            set(gca, 'LineWidth', 1.5, 'GridColor', 'black'); % Set line thickness and grid color
        end

        function plot_all_P_remnant(obj)

            pulse_time = obj.Pulse_time;
            p_remnant = obj.P_remnant;
            fitted_pt = obj.Fitted_pt;
            fitted_P = obj.Fitted_P;
            semilogx(pulse_time, p_remnant, 'o', 'MarkerSize', 8, 'Color', 'black');
            hold on;
            semilogx(fitted_pt, fitted_P, 'LineWidth',1,'LineStyle','-', 'Color','black');
            fultitle = [obj.Label_device '_' obj.Label_voltage ' fitted data'];
            title(fultitle,'FontSize', 16,'Color', 'black',Interpreter='none');
            yline(0, 'LineWidth', 1.5, 'Color', 'black'); % y = 0
            xline(0, 'LineWidth', 1.5, 'Color', 'black'); % x = 0
            xlabel('Pulse Duration (s)', 'FontSize', 16);
            ylabel('$Q_r$ ($\mu C/cm^2$)', 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'black');
            ylim([-60 60]);
            xlim([1e-8 1e-2]);
            legend([{'original data'},{'fitted data'}],Location="best");
            grid on;
            grid minor; % Use minor grid lines
            set(gca, 'LineWidth', 1.5, 'GridColor', 'black'); % Set line thickness and grid color
        end

    end

    methods(Static)

        function P_fraction = NLS_model(t,n,w,t0)
            % NLS modle
            % Ensure t > 0、tau > 0
            KAI = @(t,tau,n) (1 - exp(-(t./tau).^n));
            % 标准 Lorentzian PDF
            F = @(log_tau,w,t0) (w./((log_tau-log10(t0)).^2+w.^2))./pi;%严格按照洛伦兹分布来就不需要归一化参数A
            log_t1 = log10(t0);
            log_tau = linspace(log_t1 - 1e3*w, log_t1 + 1e3*w,1e5);
            tau =10 .^ log_tau;

            P_fraction = zeros(size(t));
            for i = 1:numel(t)
                integrand = KAI(t(i), tau, n) .* F(log_tau, w, t0);
                P_fraction(i) = trapz(log_tau, integrand);
            end
            % Normalization verification
            % norm = trapz(log_tau,F(log_tau, w, t0));
        end

    end

end

