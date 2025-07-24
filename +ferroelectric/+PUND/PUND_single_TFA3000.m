classdef PUND_single_TFA3000

    properties
        Sample_name = 'undefined';
        File_path;

        %%%%%%%%%%%%%%%%%%%%%%%% raw data %%%%%%%%%%%%%%%%%%%%%%%%
        Frequency = 0;
        Amplitude = 0;
        Pr_up = 0;
        Pr_down = 0;
        Vc_n = 0;
        Vc_p = 0;
        Thickness = 0;
        Area = 0;

        Raw_data_matrix = []; % t V I P * 5 for PUND

        V_positive = [];
        V_up = [];
        V_negative = [];
        V_down = [];

        I_positive = [];
        I_up = [];
        I_negative = [];
        I_down = [];

        T_positive = [];
        T_up = [];
        T_negative = [];
        T_down = [];
        %%%%%%%%%%%%%%%%%%%%%%%% raw data %%%%%%%%%%%%%%%%%%%%%%%%
        Q_p_raw;% unit C
        Q_u_raw;
        Q_n_raw;
        Q_d_raw;

        PUND_P; % unit uC/cm^2
        PUND_Q; % unit C
        PUND_V; % unit V
        PUND_E; % unit MV/cm

    end

    methods

        function obj = PUND_single_TFA3000(file_path,sample_name,area,thickness)
            if ~exist('sample_name','var')
                sample_name = 'undefined';
                disp(['sample_name is not defined, will extract sample_name from ',file_path]);
            end
            if ~exist('area','var')
                area = 1;
                disp(['area is not defined, will extract area from ',file_path]);
            end
            if ~exist('thickness','var')
                thickness = 1;
                disp(['thickness is not defined, will extract thickness from ',file_path]);
            end
            obj.File_path = file_path;
            obj.Sample_name = sample_name;
            obj.Area = area;
            obj.Thickness = thickness;

            obj = extract_raw_data(obj);
            obj = split_voltage(obj);
            obj = PUND_PV_cal(obj);
        end

        function obj = extract_raw_data(obj)
            file_path = obj.File_path;
            fileID = fopen(file_path, 'r');

            % Read each line using fgetl until the end of the file
            line = fgetl(fileID);
            while ischar(line)
                if strncmp(line, 'SampleName ', numel('SampleName'))
                    % Extract the next text after "SetupTitle"
                    splitLine = split(line, ':');
                    sample_name = strtrim(splitLine{2});
                elseif strncmp(line, 'Pund Frequency [Hz]', numel('Pund Frequency [Hz]'))
                    splitLine = split(line, ':');
                    val = str2double(splitLine{2});
                    frequency =  val;
                elseif strncmp(line, 'Pund Amplitude [V]', numel('Pund Amplitude [V]'))
                    splitLine = split(line, ':');
                    val = str2double(splitLine{2});
                    amplitude = val;
                elseif strncmp(line,'Pr+ [uC/cm2]',numel('Pr+ [uC/cm2]'))
                    splitLine = split(line, ':');
                    val = str2double(splitLine{2});
                    pr_up = val;
                elseif strncmp(line,'Pr- [uC/cm2]',numel('Pr- [uC/cm2]'))
                    splitLine = split(line, ':');
                    val = str2double(splitLine{2});
                    pr_down = val;
                elseif strncmp(line,'Vc- [V]',numel('Vc- [V]'))
                    splitLine = split(line, ':');
                    val = str2double(splitLine{2});
                    vc_n = val;
                elseif strncmp(line,'Vc+ [V]',numel('Vc+ [V]'))
                    splitLine = split(line, ':');
                    val = str2double(splitLine{2});
                    vc_p = val;
                elseif strncmp(line,'Thickness [nm]',numel('Thickness [nm]'))
                    splitLine = split(line, ':');
                    val = str2double(splitLine{2})*1e-9;
                    thickness = val;
                elseif strncmp(line,'Area [mm2]',numel('Area [mm2]'))
                    splitLine = split(line, ':');
                    val = str2double(splitLine{2});
                    area = val * 1e-2; % [cm2]
                elseif strncmp(line, 'Time [s]', numel('Time [s]'))
                    % 表头检测到，准备读取数据
                    data_matrix = [];
                    % 读取401行数据
                    while true
                        next_line = fgetl(fileID);
                        if ~ischar(next_line) || isempty(strtrim(next_line))
                            break;
                        end
                        val = str2double(strsplit(strtrim(next_line)));
                        data_matrix = [data_matrix; val];
                    end
                end
                % Read the next line
                line = fgetl(fileID);
            end
            fclose(fileID);
            if obj.Sample_name == "undefined"
                obj.Sample_name = sample_name;
            end
            if obj.Area == 1
                obj.Area = area;
            end
            if obj.Thickness == 1
                obj.Thickness = thickness;
            end
            obj.Frequency = frequency;
            obj.Amplitude = amplitude;
            obj.Pr_up = pr_up;
            obj.Pr_down = pr_down;
            obj.Vc_n = vc_n;
            obj.Vc_p = vc_p;
            obj.Raw_data_matrix = data_matrix;
        end

        function obj = split_voltage(obj)
            raw_data = obj.Raw_data_matrix;

            T_P = raw_data(:,1);
            T_U = raw_data(:,5);
            T_N = raw_data(:,9);
            T_D = raw_data(:,13);

            V_P = raw_data(:,2);
            V_U = raw_data(:,6);
            V_N = raw_data(:,10);
            V_D = raw_data(:,14);

            I_P = raw_data(:,3);
            I_U = raw_data(:,7);
            I_N = raw_data(:,11);
            I_D = raw_data(:,15);

            [~,idx] = findpeaks(V_P);
            P_end = 2 * idx(1);
            [~,idx] = findpeaks(V_U);
            U_end = 2 * idx(1);
            if P_end > U_end
                P_end = U_end;
            else
                U_end = P_end;
            end

            [~,idx] = findpeaks(-V_N);
            N_end = 2 * idx(1);
            [~,idx] = findpeaks(-V_D);
            D_end = 2 * idx(1);
            if N_end > D_end
                N_end = D_end;
            else
                D_end = N_end;
            end

            obj.V_positive = V_P(1:P_end);
            obj.V_up = V_U(1:U_end);
            obj.V_negative = V_N(1:N_end);
            obj.V_down = V_D(1:D_end);

            obj.I_positive = I_P(1:P_end);
            obj.I_up = I_U(1:U_end);
            obj.I_negative = I_N(1:N_end);
            obj.I_down = I_D(1:D_end);

            obj.T_positive = T_P(1:P_end);
            obj.T_up = T_U(1:U_end);
            obj.T_negative = T_N(1:N_end);
            obj.T_down = T_D(1:D_end);
        end

        function plot_PV_raw(obj,color,linewidth)
            Q_p = obj.Q_p_raw;% unit C
            Q_u = obj.Q_u_raw;
            Q_n = obj.Q_n_raw;
            Q_d = obj.Q_d_raw;

            V_p = obj.V_positive;
            V_u = obj.V_up;
            V_n = obj.V_negative;
            V_d = obj.V_down;

            area = obj.Area;% unit cm^2
            name = obj.Sample_name;

            Q_all = [Q_p;Q_u;Q_n;Q_d];
            V_all = [V_p;V_u;V_n;V_d];
            P_all = Q_all./area * 1e6;

            if ~exist("color","var")
                color = 'black';
            end
            if ~exist("linewidth","var")
                linewidth = 1.5;
            end
            plot(V_all, P_all, '-','Color',color, 'LineWidth', linewidth);
            P_max = 1.2 * max(abs(P_all));
            V_max = 1.2 * max(abs(V_all));
            xlim([-V_max,V_max]);
            ylim([-P_max,P_max]);
            title([name,'_raw_PV loop'],'FontSize',18,Interpreter='none');
            xlabel('Voltage ($V$)', 'FontSize', 12,Interpreter='latex');
            ylabel('Polarization ($\mu C/cm^2$)', 'FontSize', 12,Interpreter='latex');
        end

        function obj = PUND_PV_cal(obj)
            V_p = obj.V_positive;
            % V_u = obj.V_up;
            V_n = obj.V_negative;
            % V_d = obj.V_down;

            I_p = obj.I_positive;
            I_u = obj.I_up;
            I_n = obj.I_negative;
            I_d = obj.I_down;

            T_p = obj.T_positive;
            T_u = obj.T_up;
            T_n = obj.T_negative;
            T_d = obj.T_down;
            % T_p = obj.T_positive;
            % T_u = obj.T_positive;
            % T_n = obj.T_positive;
            % T_d = obj.T_positive;

            area = obj.Area;
            thickness = obj.Thickness;

            Q_p = cumtrapz(T_p,I_p);
            Q_u = cumtrapz(T_u,I_u);
            Q_n = cumtrapz(T_n,I_n);
            Q_d = cumtrapz(T_d,I_d);

            Q_u = Q_u+Q_p(end);
            Q_n = Q_n+Q_u(end);
            Q_d = Q_d+Q_n(end);

            Q_all = [Q_p;Q_u;Q_n;Q_d];
            Q_shift = -0.5 * (max(Q_all) + min(Q_all));

            Q_p = Q_p + Q_shift;
            Q_u = Q_u + Q_shift;
            Q_n = Q_n + Q_shift;
            Q_d = Q_d + Q_shift;
            % plot(V_p,Q_p);
            % hold on;
            % plot(V_u,Q_u);
            % hold on;
            % plot(V_n,Q_n);
            % hold on;
            % plot(V_d,Q_d);
            obj.Q_p_raw = Q_p;% unit C
            obj.Q_u_raw = Q_u;
            obj.Q_n_raw = Q_n;
            obj.Q_d_raw = Q_d;

            % cal_PUND
            I_pu = I_p-I_u;
            I_nd = I_n-I_d;

            Q_pu = cumtrapz(T_p,I_pu); % unit C
            Q_nd = cumtrapz(T_n,I_nd);
            Q_nd = Q_nd + Q_pu(end);
            Q = [Q_pu;Q_nd];
            Q = Q - 0.5*(max(Q)+min(Q));
            P = Q./area * 1e6; % unit uC/cm^2
            V = [V_p;V_n];
            E = V./thickness * 10;% unit MV/cm

            obj.PUND_Q = Q;
            obj.PUND_P = P;
            obj.PUND_V = V;
            obj.PUND_E = E;

            % plot(E,P);
        end
        function plot_PV(obj,color,linewidth,show_title)
            P = obj.PUND_P; % unit uC/cm^2
            V = obj.PUND_V; % unit V
            name = obj.Sample_name;
            if ~exist("color","var")
                color = 'black';
            end
            if ~exist("linewidth","var")
                linewidth = 1.5;
            end
            if ~exist("show_title","var")
                show_title = true;
            end
            plot(V, P, '-','Color',color, 'LineWidth', linewidth);
            P_max = 1.2 * max(abs(P));
            V_max = 1.2 * max(abs(V));
            xlim([-V_max,V_max]);
            ylim([-P_max,P_max]);
            if show_title == true
                title([name,'_PV loop'],'FontSize',18,Interpreter='none');
            end
            xlabel('Voltage ($V$)', 'FontSize', 12,Interpreter='latex');
            ylabel('Polarization ($\mu C/cm^2$)', 'FontSize', 12,Interpreter='latex');
        end


    end


end