classdef PUND_single_TFA3000

    properties
        Sample_name = '';
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

        Raw_data_matrix = [];
        %%%%%%%%%%%%%%%%%%%%%%%% raw data %%%%%%%%%%%%%%%%%%%%%%%%
        Pund_PV = [];

    end

    methods
        
        function obj = PUND_single_TFA3000(file_path)
            
            obj.File_path = file_path;
            obj = extract_raw_data(obj);
            obj = pund_PV_cal(obj,true);

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
                    area = val;
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
            obj.Sample_name = sample_name;
            obj.Frequency = frequency;
            obj.Amplitude = amplitude;
            obj.Pr_up = pr_up;
            obj.Pr_down = pr_down;
            obj.Vc_n = vc_n;
            obj.Vc_p = vc_p;
            obj.Thickness = thickness;
            obj.Area = area;
            obj.Raw_data_matrix = data_matrix;
        end

        function obj = pund_PV_cal(obj,USE_filter)
            raw_data = obj.Raw_data_matrix;
            area = obj.Area;

            T_p = raw_data(:,1);
            V_p = raw_data(:,2);
            I_p = raw_data(:,3);

            T_u = raw_data(:,5);
            V_u = raw_data(:,6);
            I_u = raw_data(:,7);

            T_n = raw_data(:,9);
            V_n = raw_data(:,10);
            I_n = raw_data(:,11);

            T_d = raw_data(:,13);
            V_d = raw_data(:,14);
            I_d = raw_data(:,15);

            if USE_filter == true
                I_p = medfilt1(I_p, 5);
                I_u = medfilt1(I_u, 5);
                I_n = medfilt1(I_n, 5);
                I_d = medfilt1(I_d, 5);
            end

            I_PU = I_p-I_u;
            I_ND = I_n-I_d;

            [V_p,I_PU,T_p] = data_clear(V_p,I_PU,T_p);
            [V_n,I_ND,T_n] = data_clear(V_n,I_ND,T_n);

            data_length = min([length(I_PU),length(I_ND)]);
            % T_n =T_p;
            % T_N = [T_p(1:data_length);T_p(data_length)+T_n(1:data_length)];
            T_N = [T_p;T_n];
            I_N = [I_PU;I_ND];

            Q = cumtrapz(T_N, I_N);

            dq = Q(data_length+1) - Q(data_length);
            Q(data_length+1:end) = Q(data_length+1:end)-dq;
            
            % Q(data_length) = ()
            Q = Q-0.5*(max(Q)+min(Q));
            voltage = [V_p;V_n];
            polarization = Q/area*1e8; % uC/cm^2
            obj.Pund_PV = [polarization voltage ];
        end

        function plot_PV(obj,color)
            P = obj.Pund_PV(:,1);
            V = obj.Pund_PV(:,2);
            plot(V, P, '-','Color',color, 'LineWidth', 1);
            xlabel('Voltage (V)', 'FontSize', 12, 'Color', 'black');
            ylabel('Polarization (\muC/cm^2)', 'FontSize', 12, 'Color', 'black');
        end
    end

    
end