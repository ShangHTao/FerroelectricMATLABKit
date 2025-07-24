classdef PUND_single_B1500
    properties
        Sample_name = '';
        CVS_file_path = '';
        Setup_title = '';
        Area = 0; % unit cm^2
        Thickness = 0; % unit nm

        %%%%%%%%%%%%%%%%%%%%%%%% raw data %%%%%%%%%%%%%%%%%%%%%%%%
        V_raw = [];
        I_raw = [];
        T_raw = [];

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

        function obj = PUND_single_B1500(cvs_file_path,sample_name,area,thickness)

            obj.CVS_file_path = cvs_file_path;
            obj.Sample_name = sample_name;
            obj.Area = area;
            obj.Thickness = thickness;

            obj = extract_raw_data(obj);
            obj = split_voltage(obj);
            obj = PUND_PV_cal(obj);
        end

        function obj = extract_raw_data(obj)
            cvs_file_path = obj.CVS_file_path;
            fileID = fopen(cvs_file_path, 'r');
            V = [];
            I = [];
            T = [];
            line = fgets(fileID);
            while ischar(line)
                if startsWith(line, 'SetupTitle, WGFMU Pattern Editor Child DataDisplay')
                    break;
                elseif startsWith(line, 'SetupTitle')
                    splitLine = split(line, ',');
                    if numel(splitLine) >= 2
                        setup_title = strtrim(splitLine{2});
                    end
                elseif startsWith(line, 'DataValue')
                    splitLine = split(line, ',');
                    V(end + 1) = str2double(splitLine{2});
                    I(end + 1) = str2double(splitLine{3});
                    T(end + 1) = str2double(splitLine{4});
                end
                line = fgetl(fileID);
            end
            fclose(fileID);
            obj.Setup_title = setup_title;
            obj.V_raw = V';
            obj.I_raw = -1*I';
            obj.T_raw = T';
        end

        function obj = split_voltage(obj,shift_n)
            V = obj.V_raw;
            I = obj.I_raw;
            T = obj.T_raw;
            VIT = [V,I,T];

            dV = diff(V);
            ddV = diff(dV);
            max_ddV = max(ddV);
            % Calculate the inflection point
            [~,idx_p] = findpeaks(ddV,'MinPeakHeight', 0.1*max_ddV);
            [~,idx_n] = findpeaks(-ddV,'MinPeakHeight', 0.1*max_ddV);

            if ~exist("shift_n","var")
                shift_n = 1;
            end

            p_start = idx_p(2)-shift_n;
            p_end = idx_p(3)+shift_n;
            u_start = idx_p(4)-shift_n;

            % Brute-force solving the optimal coordinates for u_start
            min_error = inf;
            best_u_start = u_start;
            for candidate = u_start-10 : u_start+10
                u_range = candidate : candidate + (p_end - p_start);
                if u_range(end) > size(V,1) || u_range(1) < 1
                    continue;
                end
                diff_sq = (V(p_start:p_end,:) - V(u_range,:)).^2;
                curr_error = sum(diff_sq, 'all');
                if curr_error < min_error
                    min_error = curr_error;
                    best_u_start = candidate;
                end
            end
            u_start = best_u_start;
            u_end = u_start+(p_end-p_start);% Ensure the PU length remains constant

            n_start = idx_n(5)-shift_n;
            n_end = idx_n(6)+shift_n;
            d_start = idx_n(7)-shift_n;

            % Brute-force solving the optimal coordinates for d_start
            min_error = inf;
            best_d_start = d_start;
            for candidate = d_start-10 : d_start+10
                d_range = candidate : candidate + (n_end - n_start);
                if d_range(end) > size(V,1) || d_range(1) < 1
                    continue;
                end
                diff_sq = (V(n_start:n_end,:) - V(d_range,:)).^2;
                curr_error = sum(diff_sq, 'all');
                if curr_error < min_error
                    min_error = curr_error;
                    best_d_start = candidate;
                end
            end
            d_start = best_d_start;
            d_end = d_start+(n_end-n_start); % Ensure the ND length remains constant

            VIT_p = VIT(p_start:p_end,:);
            VIT_u = VIT(u_start:u_end,:);
            VIT_n = VIT(n_start:n_end,:);
            VIT_d = VIT(d_start:d_end,:);

            obj.V_positive = VIT_p(:,1);
            obj.V_up = VIT_u(:,1);
            obj.V_negative = VIT_n(:,1);
            obj.V_down = VIT_d(:,1);

            obj.I_positive = VIT_p(:,2);
            obj.I_up = VIT_u(:,2);
            obj.I_negative = VIT_n(:,2);
            obj.I_down = VIT_d(:,2);

            obj.T_positive = VIT_p(:,3);
            obj.T_up = VIT_u(:,3);
            obj.T_negative = VIT_n(:,3);
            obj.T_down = VIT_n(:,3);
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

        function plot_IV_raw(obj,color,linewidth)
            name = obj.Sample_name;
            V = obj.V_raw;
            I = obj.I_raw;
            if ~exist("color","var")
                color = 'black';
            end
            if ~exist("linewidth","var")
                linewidth = 1.5;
            end

            plot(V, I, '-','Color',color, 'LineWidth', linewidth);
            I_max = 1.2 * max(abs(I));
            V_max = 1.2 * max(abs(V));
            xlim([-V_max,V_max]);
            ylim([-I_max,I_max]);
            title([name,'_raw_IV loop'],'FontSize',18,Interpreter='none');
            xlabel('Voltage ($V$)', 'FontSize', 12,Interpreter='latex');
            ylabel('Current ($A$)', 'FontSize', 12,Interpreter='latex');
        end
        function [Pr_p,Pr_n,Ec_p,Ec_n] = get_Pr_Ec(obj)
            P = obj.PUND_P;
            E = obj.PUND_E;

            idx_start = find(E>0);
            idx_end = find(E<0);
            E = E(idx_start(1):idx_end(end));
            P = P(idx_start(1):idx_end(end));

            [x_c,y_c] = find_axis_crossings(E, P);

            idx = x_c > 0;
            Ec_p = mean(x_c(idx));% get mean

            idx = x_c < 0;
            Ec_n = mean(x_c(idx));

            idx = y_c > 0;
            Pr_p = mean(y_c(idx));

            Pr_n = 0.5 * (P(1) + P(2));
        end

    end


end