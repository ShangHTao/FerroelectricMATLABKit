classdef switching_dynamic_multi_voltage
    % 多电压SD测试
    properties
        N_voltage;
        Voltage_values = [];

        Thickness = 0;% unit nm
        Label_voltage = {};
        Label_device = '';

        SV_objs = {};% SV - > single_voltage -> single_voltage_objs

        Fitted_t0 = 0;
        Fitted_w = 0;
        KB_t0 = 0;
        KB_w = 0;
        Folder_path = '';
    end

    methods
        function obj = switching_dynamic_multi_voltage(folder_path,area,thickness,v_bias)
            parts = split(folder_path, filesep);
            contents = dir(folder_path);
            isSubFolder = [contents.isdir] & ~ismember({contents.name}, {'.', '..'}) & contains({contents.name}, 'V');
            subFolders = contents(isSubFolder);
            voltageValues = zeros(1, numel(subFolders));% 提取每个文件夹名中的电压数值部分
            for i = 1:numel(subFolders)
                name = subFolders(i).name;
                voltageValues(i) = sscanf(name, '%f');
            end
            [~, sortIdx] = sort(voltageValues);% 按照电压值排序
            subFolders = subFolders(sortIdx);
            N = numel(subFolders);
            objs = cell(N,1);
            label_voltage = cell(N,1);
            for i = 1:N

                SV_folder_path = [subFolders(i).folder filesep subFolders(i).name];
                objs{i} = switching_dynamic_single_voltage(SV_folder_path,area,thickness,v_bias);
                label_voltage{i} = subFolders(i).name;
            end
            obj.Label_device = parts(end-1);
            obj.Folder_path = folder_path;
            obj.Voltage_values = voltageValues(sortIdx)';
            obj.Label_voltage = label_voltage;
            obj.N_voltage = N;
            obj.SV_objs = objs;
            obj.Thickness = thickness;

            obj = plot_Merz_cal(obj,false);
        end

        function plot_P_remnant_mv(obj)
            % % figure;
            N = obj.N_voltage;
            color = hsv(N);
            legend_str = cell(1,N);
            for i = 1:N
                sv_obj = obj.SV_objs{i};
                legend_str{i} =  [sv_obj.Label_device '_' sv_obj.Label_voltage];
                pulse_time = sv_obj.Pulse_time;
                p_remnant = sv_obj.P_remnant;
                semilogx(pulse_time, p_remnant, 'o', 'MarkerSize', 8, 'Color', color(i,:));
                hold on;
            end
            fultitle = obj.Label_device;
            title(fultitle,'FontSize', 16,'Color', 'black',Interpreter='none');
            yline(0, 'LineWidth', 1.5, 'Color', 'black'); % y = 0
            xline(0, 'LineWidth', 1.5, 'Color', 'black'); % x = 0
            xlabel('Pulse Duration (s)', 'FontSize', 16);
            ylabel('$Q_r$ ($\mu C/cm^2$)', 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'black');
            ylim([-60 60]);
            xlim([1e-8 1e-2]);
            legend(legend_str,'Location','best',Interpreter='none');
            grid on;
            grid minor; % Use minor grid lines
            set(gca, 'LineWidth', 1.5, 'GridColor', 'black'); % Set line thickness and grid color
        end

        function plot_fitted_P_remnant_mv(obj)
            N = obj.N_voltage;
            color = hsv(N);
            legend_str = cell(1,N);
            for i = 1:N
                sv_obj = obj.SV_objs{i};
                legend_str{i} =  [sv_obj.Label_device '_' sv_obj.Label_voltage];
                fitted_logpt = sv_obj.Fitted_pt;
                fitted_P = sv_obj.Fitted_P;
                semilogx(fitted_logpt, fitted_P, 'LineWidth',1,'LineStyle','-', 'Color', color(i,:));
                hold on;
            end
            fultitle = [obj.Label_device ' fitted data'];
            title(fultitle,'FontSize', 16,'Color', 'black',Interpreter='none');
            yline(0, 'LineWidth', 1.5, 'Color', 'black'); % y = 0
            xline(0, 'LineWidth', 1.5, 'Color', 'black'); % x = 0
            xlabel('Pulse Duration (s)', 'FontSize', 16);
            ylabel('$Q_r$ ($\mu C/cm^2$)', 'Interpreter', 'latex', 'FontSize', 16, 'Color', 'black');
            ylim([-60 60]);
            xlim([1e-8 1e-2]);
            legend(legend_str,'Location','best',Interpreter='none');
            grid on;
            grid minor; % Use minor grid lines
            set(gca, 'LineWidth', 1.5, 'GridColor', 'black'); % Set line thickness and grid color
        end

        function plot_all_P_remnant_mv(obj)
            N = obj.N_voltage;
            color = hsv(N);
            legend_str = cell(1,2*N);
            for i = 1:N
                sv_obj = obj.SV_objs{i};
                legend_str{i} =  [sv_obj.Label_device '_' sv_obj.Label_voltage];
                pulse_time = sv_obj.Pulse_time;
                p_remnant = sv_obj.P_remnant;
                semilogx(pulse_time, p_remnant, 'o', 'MarkerSize', 8, 'Color', color(i,:));
                hold on;
            end
            for i = 1:N
                sv_obj = obj.SV_objs{i};
                legend_str{N+i} =  [sv_obj.Label_device '_' sv_obj.Label_voltage '_fitted'];
                fitted_logpt = sv_obj.Fitted_pt;
                fitted_P = sv_obj.Fitted_P;
                semilogx(fitted_logpt, fitted_P, 'LineWidth',1,'LineStyle','-', 'Color', color(i,:));
                hold on;
            end
            fultitle = [obj.Label_device];
            title(fultitle,'FontSize', 16,'Color', 'black',Interpreter='none');
            yline(0, 'LineWidth', 1.5, 'Color', 'black'); % y = 0
            xline(0, 'LineWidth', 1.5, 'Color', 'black'); % x = 0
            xlabel('Pulse Duration (s)', 'FontSize', 16,Interpreter='latex');
            ylabel('$Q_r$ ($\mu C/cm^2$)','FontSize', 16, 'Color', 'black',Interpreter='latex');
            ylim([-60 60]);
            xlim([1e-8 1e-2]);
            legend(legend_str,'Location','best',Interpreter='none');
            grid on;
            grid minor; % Use minor grid lines
            set(gca, 'LineWidth', 1.5, 'GridColor', 'black'); % Set line thickness and grid color
        end

        function plot_PVT(obj)

            N_v = obj.N_voltage;
            N_pt = 0;
            for i = 1:N_v %求出脉冲时间点数
                N = obj.SV_objs{i}.N_test;
                if N_pt < N
                    N_pt = N;
                    logt = log10(obj.SV_objs{i}.Pulse_time);
                end
            end
            %脉冲线性插值(由于1V存在两个坏点) 2*N_pt
            logt_max = max(logt);
            logt_min = min(logt);
            u_logt = linspace(logt_min,logt_max,2*N_pt);
            for i = 1:N_v
                p_i = obj.SV_objs{i}.P_remnant';
                logt_i = log10(obj.SV_objs{i}.Pulse_time)';
                p_inter(i,:) = interp1(logt_i, p_i, u_logt, 'linear', 'extrap');
            end
            voltage = obj.Voltage_values;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % imagesc(u_logt,voltage,p_inter);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            imagesc(u_logt,voltage,abs(p_inter));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            title(obj.Label_device,'FontSize',16);
            ylabel('applied voltage (V)', 'FontSize', 16,Interpreter='latex');
            xlabel('$log_{10} (t) \; (log_{10}(s))$', 'FontSize',16,Interpreter='latex');
            c = colorbar();
            c.Label.String = 'remnant polarization $(\mu C / cm^2)$';
            c.Label.Interpreter  = 'latex';
            c.Label.FontSize = 16;
        end

        function obj = plot_Merz_cal(obj,draw)
            % figure;
            N_v = obj.N_voltage;
            fitted_t0 = zeros(1,N_v);
            fitted_w = zeros(1,N_v);
            th = zeros(1,N_v);
            for i = 1:N_v
                sv_obj = obj.SV_objs{i};
                fitted_t0(i) = sv_obj.Fitted_t0;
                fitted_w(i) = sv_obj.Fitted_w;
                th(i) = sv_obj.Thickness;
            end
            logt0 = log10(fitted_t0);

            voltage = obj.Voltage_values';
            E_b = 1./(voltage./th); %计算电场的倒数
            E_b2 = E_b.^2;
            % 最小二乘法求比例系数k = (x * y') / (x * x');
            t0_k_b = polyfit(E_b, logt0, 1);
            w0_k_b = polyfit(E_b2, fitted_w, 1);
            % k_t0 = (E_b*fitted_t0')/(E_b*E_b');
            % k_w = (E_b2*fitted_w')/(E_b2*E_b2');

            f_E_b = linspace(min(E_b),max(E_b),100);
            f_logt0 = t0_k_b(1) .* f_E_b + t0_k_b(2);

            f_E_b2 = linspace(min(E_b2),max(E_b2),100);
            f_w = w0_k_b(1) .* f_E_b2 + w0_k_b(2);

            obj.Fitted_t0 = fitted_t0;
            obj.Fitted_w = fitted_w;
            obj.KB_t0 = t0_k_b;
            obj.KB_w = w0_k_b;

            if draw == true
                subplot(1,2,1);
                plot(E_b,logt0,'o', 'MarkerSize', 8, 'Color','r');
                hold on;
                plot(f_E_b,f_logt0,'--', 'LineWidth', 2, 'Color','g');
                xlim([min(E_b) max(E_b)]);
                ylim([min(logt0) max(logt0)]);
                xlabel('1/E (nm V^{-1})','FontSize',16,Interpreter='tex');
                ylabel('t_{mean} (s)','FontSize',16,Interpreter='tex');
                title('mean','FontSize',16,Interpreter='latex');

                subplot(1,2,2);
                plot(E_b2,fitted_w,'o', 'MarkerSize', 8, 'Color','r');
                hold on;
                plot(f_E_b2,f_w,'--', 'LineWidth', 2, 'Color','g');
                xlim([min(E_b2) max(E_b2)]);
                ylim([min(fitted_w) max(fitted_w)]);
                xlabel('1/E^2 (nm^2 V^{-2})','FontSize',16,Interpreter='tex');
                ylabel('w (s)','FontSize',16,Interpreter='tex');
                title('HWHM','FontSize',16,Interpreter='tex');
            end
        end

        function plot_PVT_fitted(obj)
            N_v = obj.N_voltage;
            N_int_voltage = 100;

            th = obj.Thickness;

            voltage = obj.Voltage_values;
            voltage_start = voltage(1);
            voltage_stop = voltage(end);
            voltage_int = linspace(voltage_start,voltage_stop,N_int_voltage);% 生成100个电压电

            P_max = zeros(1,N_v);

            for i = 1:N_v
                sv_obj = obj.SV_objs{i};
                polarity = sv_obj.Polarity(1);
                P_max(i) = max(abs(sv_obj.P_remnant)) * polarity; %电压正负号也加入
            end
            P_max_int = interp1(voltage,P_max,voltage_int,"linear","extrap");


            fitted_pt = obj.SV_objs{1}.Fitted_pt;
            kb_t0 = obj.KB_t0;
            kb_w = obj.KB_w;

            E_b = 1./(voltage_int./th);
            E_b2 = E_b.^2;
            logt0 = kb_t0(1) .* E_b + kb_t0(2);
            w = kb_w(1) .* E_b2 + kb_w(2);

            fraction_P_fit = zeros(N_int_voltage,numel(fitted_pt));
            P_int = zeros(N_int_voltage,numel(fitted_pt));
            for i = 1:numel(voltage_int)
                fraction_P_fit(i,:) = switching_dynamic_single_voltage.NLS_model(fitted_pt,2, w(i), 10^logt0(i));
                P_int(i,:) = fraction_P_fit(i,:) .* P_max_int(i);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % imagesc(log10(fitted_pt),voltage_int,P_int)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            imagesc(log10(fitted_pt),voltage_int,abs(P_int));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            title(obj.Label_device,'FontSize',16);
            ylabel('applied voltage (V)', 'FontSize', 16,Interpreter='latex');
            xlabel('$log_{10} (t) \; (log_{10}(s))$', 'FontSize',16,Interpreter='latex');
            c = colorbar();
            c.Label.String = 'remnant polarization $(\mu C / cm^2)$';
            c.Label.Interpreter  = 'latex';
            c.Label.FontSize = 16;
        end
    end
end