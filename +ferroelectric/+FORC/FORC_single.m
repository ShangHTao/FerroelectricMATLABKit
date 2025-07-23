classdef FORC_single
    properties
        CVS_file_path;
        File_name;
        Sample_name = 'undefined';
        Apply_voltage; % unit V
        Thickness = 0; % unit nm
        Area = 0; % unit cm^2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% raw data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Setup_title = '';
        Voltage_raw = [];
        Time_raw = [];
        Current_raw = [];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% raw data %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% processed data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_triangular_wave = 0;
        Voltage_reversal = [];
        Voltage = [];
        Time = [];
        Current = [];

        Voltage_split = {}; %三角波电压划分
        Time_split = {};
        Current_split = {};
        Polarization_split = {}; % unit uC/cm^2

        %存储为矩阵更为合理
        Voltage_split_reversal = []; %三角波回转电压划分
        Time_split_reversal = [];
        Current_split_reversal = [];
        Polarization_split_reversal = []; % unit uC/cm^2

        Rho_raw_v = []; % unit uC/(MV*nm)^2
        Rho_raw_E = []; % unit uC/MV^2
        V_grid_raw = [];% unit V
        Vr_grid_raw = [];% unit V


        % FORC插值后的数据
        Rho_E_q = []; % unit uC/(MV*nm)^2
        V_grid_q = []; % unit V
        Vr_grid_q = []; % unit V

        Rho_E_trans = []; %  Coordinate transformation
        Vc_grid = [];
        Vb_grid = [];

        Const_dV_dt = 0;
        Const_dt = 0;
        Const_dV = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% processed data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    properties (Access=private)

        split_voltage_th = 0.9;
        current_compensation = true;
        cutEdgePoint = 1;
        cutRampPoints = 1;
    end

    methods
        function obj = FORC_single(cvs_file_path,area,thickness,apply_voltage,args)
            arguments
                cvs_file_path (1, :) char
                area (1,1) double
                thickness (1,1) double
                apply_voltage (1,1) double
            end

            arguments
                args.Sample_name = "undefined";
                args.split_voltage_th = 0.9;
                args.current_compensation = true;
                args.cutEdgePoint = 1;
                args.cutRampPoints = 1;
            end

            if exist(cvs_file_path,"file") ~= 2
                error(['There is no file named',cvs_file_path ]);
            end
            [~, name, ext] = fileparts(cvs_file_path);
            filename = [name, ext];

            obj.CVS_file_path = cvs_file_path;
            obj.Area = area;
            obj.Thickness = thickness;
            obj.Apply_voltage = apply_voltage;
            obj.File_name = filename;
            obj.Sample_name = args.Sample_name;
            obj.split_voltage_th = args.split_voltage_th;
            obj.current_compensation = args.current_compensation;
            obj.cutEdgePoint = args.cutEdgePoint;
            obj.cutRampPoints = args.cutRampPoints;

            obj = extract_raw_data(obj);
            obj = split_voltage(obj,obj.split_voltage_th);
            obj = cal_PV_loop(obj,obj.current_compensation);
            obj = split_reversal_voltage(obj,obj.cutEdgePoint);
            obj = cal_FORC_use_current(obj,obj.cutRampPoints);
            obj = interpolate_FORC(obj);
            obj = coordinate_trans(obj);

        end

        function obj = extract_raw_data(obj)
            cvs_file_path = obj.CVS_file_path;
            fileID = fopen(cvs_file_path, 'r');
            % 第一遍读取确定有效行数据
            nLines = 0;
            while ~feof(fileID)
                line = fgetl(fileID);
                if startsWith(line, 'DataValue')
                    nLines = nLines + 1;
                elseif startsWith(line, 'SetupTitle, WGFMU Pattern Editor Child DataDisplay')
                    break; % 结束读取后面重复数据
                end
            end
            fclose(fileID);

            V = zeros(nLines, 1);
            I = zeros(nLines, 1);
            T = zeros(nLines, 1);
            % 第二遍读取，提取数据
            fileID = fopen(cvs_file_path, 'r');
            index = 0;
            while ~feof(fileID)
                line = fgetl(fileID);
                if startsWith(line, 'SetupTitle, WGFMU Pattern Editor Child DataDisplay')
                    break;  % 结束读取后面重复数据
                elseif startsWith(line, 'SetupTitle')
                    splitLine = split(line, ',');
                    if numel(splitLine) >= 2
                        setup_title = strtrim(splitLine{2});
                    end
                elseif startsWith(line, 'DataValue')
                    splitLine = split(line, ',');
                    index = index + 1;
                    V(index) = str2double(splitLine{2});
                    I(index) = str2double(splitLine{3});
                    T(index) = str2double(splitLine{4});
                    % raw_data{end+1} = splitLine;  % 存储为 cell 的 cell，每行一个单元
                end
            end
            fclose(fileID);

            obj.Setup_title = setup_title;
            obj.Voltage_raw = V;
            obj.Current_raw = -1*I;
            obj.Time_raw = T;
        end

        function plot_raw_VIT(obj)
            set_title = obj.Setup_title;
            t_us = obj.Time_raw/1e-6;
            v = obj.Voltage_raw;
            c_uA = obj.Current_raw/1e-6;
            yyaxis left;
            plot(t_us, v,'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'none', 'MarkerSize', 8, 'Color', 'blue');
            % Set appropriate labels and title for the left y-axis plot
            xlabel('Time (\mus)', 'FontSize',12);
            ylabel('Voltage (V)', 'FontSize',12, 'Color', 'blue');
            set(gca, 'FontSize', 12);  % Set font size
            % Automatically adjust the y-axis to be symmetric around zero
            ylimit = max(abs(ylim));
            ylim([-ylimit ylimit]);
            hold on; % Enable holding the plot
            yyaxis right;
            plot(t_us, c_uA,'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 4, 'Color', 'red');
            % Set appropriate label for the right y-axis plot
            ylabel('Current (\muA)', 'FontSize',12,'Color', 'red');
            %-----------------raw TVI curve-----------------%
            grid on;
            grid minor;  % Use minor grid lines
            set(gca, 'LineWidth', 1.5, 'GridColor', 'black');  % Set line thickness and color
            % Set the legend for both plots
            legend('V-t', 'I-t', 'FontSize',12,'Location','Best');
            title(['VIT--',set_title],'FontSize', 12,Interpreter='none');
        end

        function plot_raw_VI(obj)
            set_title = obj.Setup_title;
            v = obj.Voltage_raw;
            c_uA = obj.Current_raw/1e-6;
            plot(v, c_uA,'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'none', 'MarkerSize', 8, 'Color', [0 0.4470 0.7410]);
            xline(0,'LineWidth', 1.5,'Color', 'black'); % x = 0
            yline(0,'LineWidth', 1.5,'Color', 'black'); % y = 0
            xlabel('Voltage (V)', 'FontSize',12);
            ylabel('Current (\muA)', 'FontSize',12,'Color', 'black');
            legend('I-V', 'FontSize',12, 'Location','Best');
            ylimit = max(abs(ylim));
            ylim([-ylimit ylimit]);
            xlimit = max(abs(xlim));
            xlim([-xlimit xlimit])
            set(gca, 'FontSize', 12);  % Set font size
            grid on;
            grid minor;  % Use minor grid lines
            set(gca, 'LineWidth', 1.5, 'GridColor', 'black');  % Set line thickness and color
            title(['I-V loop--',set_title],'FontSize', 12);
        end

        function obj = split_voltage(obj,split_voltage_th)
            % split_voltage_th Separation voltage threshold coefficient. default split_voltage_th = 0.9
            % if ~exist('split_voltage_th', 'var') || isempty(split_voltage_th)
            %     split_voltage_th = 0.9;
            % end
            v = obj.Voltage_raw;
            t = obj.Time_raw;
            c = obj.Current_raw;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % Exclude the first few data points
            % % Exclude the last few data points
            % dv = diff(v);
            % idx = find(dv<0);
            % n1 = idx(1);
            % dv = dv(n1:end);
            % idx = find(dv>0);
            % n2 = idx(1);
            % n3 = n1 + n2;
            % dv = dv(n2:end);
            % idx = find(dv(1:end-1).*dv(2:end) < 0); % 寻找过零点
            % n4 = n3+idx(end);
            % v = v(n3:n4);
            % t = t(n3:n4);
            % t = t-t(1);
            % c = c(n3:n4);
            dv = diff(v);
            idx = find(dv<0);
            n1 = idx(1);
            while v(n1)-v(n1+1) <= 0.01
                n1 = n1+1;
            end
            dv = dv(n1:end);
            idx = find(dv(1:end-1).*dv(2:end) < 0); % 寻找过零点
            n2 = n1+idx(end);
            v = v(n1:n2);
            t = t(n1:n2);
            c = c(n1:n2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %求回转电压Vr
            [~, locs_max] = findpeaks(v);
            [~, locs_min] = findpeaks(-v);
            Vr = v(locs_min);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %分离每个三角波电压
            %顺便求取dV/dt
            dv_max = max(abs(dv));
            dv_satisfy = abs(dv)>split_voltage_th*dv_max;
            dv_const = mean(abs(dv(dv_satisfy)));
            dt_const = mean(diff(t));
            % fprintf('const dv/dt = %f V/s\n',dv_const/dt_const);
            N = numel(Vr);
            sp_v = cell(N,1);
            sp_c = cell(N,1);
            sp_t = cell(N,1);
            locs_max = [1;locs_max;numel(v)];%Add 1 at the start of the list
            for i = 1:N
                sp_v{i} = v(locs_max(i):locs_max(i+1));
                sp_c{i} = c(locs_max(i):locs_max(i+1));
                sp_t{i} = t(locs_max(i):locs_max(i+1));
                % %%%
                % figure;
                % plot(sp_t{i}/1e-6,sp_v{i},'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'none', 'MarkerSize', 8, 'Color', 'blue');
                % % Set appropriate labels and title for the left y-axis plot
                % xlabel('Time (\mus)', 'FontSize',12);
                % ylabel('Voltage (V)', 'FontSize',12, 'Color', 'blue');
                % set(gca, 'FontSize', 12);  % Set font size
                % % Automatically adjust the y-axis to be symmetric around zero
                % ylimit = max(abs(ylim));
                % ylim([-ylimit ylimit]);
                % hold on; % Enable holding the plot
                % yyaxis right;
                % plot(sp_t{i}/1e-6,sp_c{i}/1e-6,'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 4, 'Color', 'red');
                % % Set appropriate label for the right y-axis plot
                % ylabel('Current (\muA)', 'FontSize',12,'Color', 'red');
                % grid on;
                % grid minor;  % Use minor grid lines
                % set(gca, 'LineWidth', 1.5, 'GridColor', 'black');  % Set line thickness and color
                % % Set the legend for both plots
                % legend('V-t', 'I-t', 'FontSize',12,'Location','Best');
                % %%%
            end
            obj.Voltage = v;
            obj.Current = c;
            obj.Time = t;
            obj.Voltage_split = sp_v;
            obj.Current_split = sp_c;
            obj.Time_split = sp_t;
            obj.N_triangular_wave = N;
            obj.Voltage_reversal = Vr;
            obj.Const_dV = dv_const;
            obj.Const_dt = dt_const;
            obj.Const_dV_dt = dv_const/dt_const;
        end

        function obj = cal_PV_loop(obj,current_compensation)
            % v_sp = obj.Voltage_split;
            c_sp = obj.Current_split; % unit A
            t_sp = obj.Time_split; % unit s
            N = obj.N_triangular_wave;
            area = obj.Area; % unit cm^2

            p = cell(N,1);
            if current_compensation == true
                for i = 1:N
                    I = c_sp{i};
                    t = t_sp{i};
                    Q = cumtrapz(t, I); % unit C

                    I_bias_step = 1e-6 * max(abs(I));
                    bias_val = Q(1) - Q(end);
                    n_bias = 1;
                    condition_met = abs(Q(1) - Q(end)) / max(abs(Q)) < 0.01; % 闭口偏差小于1%
                    while ~condition_met
                        if bias_val > 0
                            Q_sign = 1;
                        else
                            Q_sign = -1;
                        end

                        I = I + Q_sign * I_bias_step * n_bias;
                        Q = cumtrapz(t,I);

                        bias_val = Q(1) - Q(end);
                        n_bias = n_bias + 1;
                        condition_met = abs(Q(1) - Q(end)) / max(abs(Q)) < 0.01; % 闭口偏差小于1%
                        if n_bias > 1e4 % 迭代超过3次就不收敛放弃电流补偿
                            Q = cumtrapz(t,c_sp{i});
                            break;
                        end
                    end
                    p{i} = Q/area * 1e6; % unit uC/cm^2
                    % plot(v_sp{i},p{i});
                    % hold on;
                end
            else
                for i = 1:N
                    I = c_sp{i};
                    t = t_sp{i};
                    Q = cumtrapz(t,I);
                    p{i} = Q/area * 1e6; % unit uC/cm^2
                    % plot(v_sp{i},p{i});
                    % hold on;
                end

            end
            obj.Polarization_split = p;
        end

        function obj = split_reversal_voltage(obj,cutEdgePoint)
            sp_v = obj.Voltage_split;
            sp_c = obj.Current_split;
            sp_t = obj.Time_split;
            sp_p = obj.Polarization_split;
            v_r = obj.Voltage_reversal;
            N = obj.N_triangular_wave;
            % if ~exist('cutEdgePoint','var')
            %     cutEdgePoint = 1; % default = 1
            % end

            idx = find(sp_v{end} == v_r(end));
            N_reversal = numel(sp_v{end}(idx:end))-cutEdgePoint; % -1 目的同下
            % N_v = numel(sp_v{end}(idx:end));

            sp_vr = zeros(N,N_reversal);
            sp_cr = zeros(N,N_reversal);
            sp_tr = zeros(N,N_reversal);
            sp_pr = zeros(N,N_reversal);

            for i = 1:N
                idx = find(sp_v{i} == v_r(i));
                sp_vr_i = sp_v{i}(idx:end-cutEdgePoint);% ！！！！！！
                sp_cr_i = sp_c{i}(idx:end-cutEdgePoint);% end-cutEdgePoint 的目的是去除电压为负值的无用信号，FORC图像边缘会出现
                sp_tr_i = sp_t{i}(idx:end-cutEdgePoint);% 通常来说cutEdgePoint = 即可
                sp_pr_i = sp_p{i}(idx:end-cutEdgePoint);

                N_sp_vr_i = numel(sp_vr_i);

                sp_vr(i,end-N_sp_vr_i+1:end) = sp_vr_i';
                sp_cr(i,end-N_sp_vr_i+1:end) = sp_cr_i';
                sp_tr(i,end-N_sp_vr_i+1:end) = sp_tr_i';
                sp_pr(i,end-N_sp_vr_i+1:end) = sp_pr_i';

                % plot(sp_vr_i,sp_cr_i);
                % plot(sp_tr_i,sp_vr_i);
                % hold on;
            end
            obj.Voltage_split_reversal = sp_vr;
            obj.Current_split_reversal = sp_cr;
            obj.Time_split_reversal = sp_tr;
            obj.Polarization_split_reversal = sp_pr;
        end

        function obj = cal_FORC_use_current(obj,cutRampPoints)
            % if ~exist('cutRampPoints','var')
            %     cutRampPoints = 1; % default = 1
            % end

            % rho(Er,E) = 0.5 * partial2(P)/ (partial(Er) * partial(E))
            % ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
            % rho(Er,E) = 0.5 * 1/dE_dt * J(Eri, E)-J(Eri_1, E) / (Eri-Eri_1)
            % E = dV/dl assume a constant electric field
            % Replace E with V
            % define rho_v(Vr,V) = 0.5 * 1/dV_dt * J(Vri, V)-J(Vri_1, V) / (Vri-Vri_1)
            % rho(Er,E) = l^2 * rho_v
            dv_dt = obj.Const_dV_dt; % unit V/s
            area = obj.Area; % unit cm^2
            c = obj.Current_split_reversal; % unit A
            j = c./area;% unit A/cm^2;
            v_grid = obj.Voltage_split_reversal;
            vr = obj.Voltage_reversal;
            thickness = obj.Thickness;
            % 补齐aplied v
            [rows, cols] = size(v_grid);
            for i = 1:cols
                idx_all = find(v_grid(:, i) == 0);
                if ~isempty(idx_all)
                    idx = max(idx_all);
                    if idx < rows
                        v_mean = mean(v_grid(idx+1:end, i), 'omitnan');
                        v_grid(1:idx, i) = v_mean;
                    else
                        % 如果0出现在最后一行，不做处理或设为默认值
                        % v_grid(1:idx, i) = some_default_value; % 如0或NaN
                    end
                end
            end
            [rows, cols] = size(j);
            djdvr = zeros(rows,cols); % unit A/(V*cm^2);

            for i = 2:rows
                j_i = j(i,:);
                j_i_1 = j(i-1,:);
                % idx = find(j_i ~= 0);
                % j_i([idx(1) idx(1)+edgeCutPoints]) = 0;
                % j_i_1(idx(1)+edgeCutPoints) = 0;

                djdvr(i,:) = (j_i-j_i_1)/(vr(i)-vr(i-1));
                idx = find(djdvr(i,:) ~= 0);

                if cutRampPoints >= numel(idx)
                    djdvr(i,:) = 0;
                else
                    djdvr(i,idx(1):idx(1)+cutRampPoints) = 0;
                end
            end

            rho_raw_v =  0.5 * (1/dv_dt) * djdvr * 1e4;  % unit uC/(MV*nm)^2
            rho_raw_E = thickness^2 * rho_raw_v; % unit uC/MV^2
            vr_grid = repmat(vr, 1, size(v_grid, 2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % surf(v_grid,vr_grid,rho_raw_v);% 可以绘制一些了
            % view(2);
            % colormap("hot");

            obj.Rho_raw_v = rho_raw_v;
            obj.Rho_raw_E = rho_raw_E;
            obj.V_grid_raw = v_grid;
            obj.Vr_grid_raw = vr_grid;
        end

        function plot_FORC_raw(obj)
            rho_raw_E = abs(obj.Rho_raw_E); % unit uC/MV^2
            thickness = obj.Thickness;
            vr_grid = obj.Vr_grid_raw ./ thickness * 1e1; % MV/cm
            v_grid = obj.V_grid_raw ./ thickness * 1e1;

            surf(v_grid,vr_grid,rho_raw_E);
            view(2);
            colormap("jet");
            xlabel('E(MV/cm^2)',LineWidth=8,Interpreter='tex');
            ylabel('E_r(MV/cm^2)',LineWidth=8,Interpreter='tex');
        end

        function obj = interpolate_FORC(obj,N_order)
            % v grid 不均匀，需要进行插值处理
            % 处理为N阶的方阵 default = 500
            if ~exist('N_order', 'var') || isempty(N_order)
                N_order = 500;
            end
            rho_E = obj.Rho_raw_E;
            v_grid = obj.V_grid_raw;
            vr_grid = obj.Vr_grid_raw;
            V_app = obj.Apply_voltage;

            [row_rho, col_rho] = size(rho_E);
            [row_v, col_v] = size(v_grid);
            [row_vr, col_vr] = size(vr_grid);

            if row_rho ~= row_v && row_v ~= row_vr && col_rho ~= col_v && col_v ~= col_vr
                error('The dimensions of the three matrices rho_raw_E, v_grid, and vr_grid are inconsistent.');
            end

            %采取两种插值方式，先对每一行进行插值
            rho_E_row_q = zeros(row_rho,N_order); % row_rho * N_order的矩阵
            v_grid_row_q = zeros(row_rho,N_order);
            for i = 1:row_rho
                v_i = linspace(-V_app,V_app,N_order);

                rho_i = interp1(v_grid(i,:),rho_E(i,:),v_i,"linear","extrap");

                rho_E_row_q(i,:) = rho_i;
                v_grid_row_q(i,:) = v_i;
            end

            % 现在对每一列进行插值,二维插值
            % 由于Vr的特殊性，对Vr拓展值
            vr_grid =repmat(vr_grid(:,1), 1, N_order); % size = 38*N
            v = linspace(-V_app,V_app,N_order);

            [v_grid_q, vr_grid_q] = meshgrid(v, -v); % size = N*N '-'保证vr符号的一致性

            % rho_E_q = zeros(N_order);
            rho_E_q = interp2(v_grid_row_q,vr_grid,rho_E_row_q,v_grid_q,vr_grid_q,"linear",0); % NaN替换为0

            obj.Rho_E_q = rho_E_q;
            obj.V_grid_q = v_grid_q;
            obj.Vr_grid_q = vr_grid_q;
        end

        function obj = coordinate_trans(obj)
            v = obj.V_grid_q;
            vr = obj.Vr_grid_q;
            rho = obj.Rho_E_q;

            N_grid = size(rho,1);

            v_flat  = v(:);
            vr_flat = vr(:);
            rho_flat = rho(:);

            vc_flat = 0.5 * (v_flat - vr_flat); % Vc -> V coercive
            vb_flat = 0.5 * (v_flat + vr_flat); % Vb -> V bias

            vcq_range = linspace(min(vc_flat), max(vc_flat), N_grid);
            vbq_range = linspace(min(vb_flat), max(vb_flat), N_grid);
            [vcq, vbq] = meshgrid(vcq_range, vbq_range);
            rho_q = griddata(vc_flat, vb_flat, rho_flat, vcq, vbq, 'linear');
            rho_q(isnan(rho_q)) = 0;

            % imagesc(vcq_range,vbq_range,rho_q);
            % colorbar;
            % xlabel('Vc'); ylabel('Vb');
            % set(gca,'YDir','normal');% Reverse the Y-axis

            obj.Rho_E_trans = rho_q;
            obj.Vc_grid = vcq;
            obj.Vb_grid = vbq;
        end

        function plot_FORC_inter(obj,color,use_electric_field,show_title)
            rho_E = abs(obj.Rho_E_q);
            v = obj.V_grid_q(1,:);% unit V
            vr = obj.Vr_grid_q(:,1);
            thickness = obj.Thickness; % unit nm

            if use_electric_field == true
                x = v./thickness * 1e1; % unit MV/cm
                y = vr./thickness * 1e1;
                imagesc(x,y,rho_E);
                xlabel('E (MV/cm)', 'FontSize', 12);
                ylabel('E_r (MV/cm)', 'FontSize', 12);
                set_title = 'FORC Distribution \rho(E, E_r)';

            else
                x = v;% unit V
                y = vr;
                imagesc(x,y,rho_E);
                xlabel('Voltage (V)', 'FontSize', 12);
                ylabel('Voltage_r (V)', 'FontSize', 12);
                set_title = 'FORC Distribution \rho(V, V_r)';
            end
            if show_title == true
                title(set_title,'FontSize', 12);
            end
            colormap(color);
            % Set equal aspect ratio to make the colormap square
            pbaspect([1 1 1]);
            % Reverse the Y-axis
            set(gca,'YDir','normal','FontSize', 12);
            % Ensure all integers within xlim and ylim are shown
            xticks(floor(min(xlim)):ceil(max(xlim)));
            yticks(floor(min(ylim)):ceil(max(ylim)));
        end

        function plot_FORC_trans(obj,color,use_electric_field,show_title,show_max)
            if ~exist('color','var')
                color = jet;
            end
            if ~exist('use_electric_field','var')
                use_electric_field = false;
            end
            if ~exist('show_title','var')
                show_title = false;
            end
            if ~exist('show_max','var')
                show_max = false;
            end

            rho_E = abs(obj.Rho_E_trans);
            vc = obj.Vc_grid(1,:);% unit V
            vb = obj.Vb_grid(:,1);
            thickness = obj.Thickness; % unit nm

            if use_electric_field == true
                x = vc./thickness * 1e1; % unit MV/cm
                y = vb./thickness * 1e1;
                imagesc(x,y,rho_E);
                xlabel('E_c (MV/cm)', 'FontSize', 12);
                ylabel('E_{bias} (MV/cm)', 'FontSize', 12);
                set_title = 'FORC Distribution \rho(E_c, E_{bias})';
            else
                x = vc;% unit V
                y = vb;
                imagesc(x,y,rho_E);
                xlabel('Voltage_c (V)', 'FontSize', 12);
                ylabel('Voltage_{bias} (V)', 'FontSize', 12);
                set_title = 'FORC Distribution \rho(V_c, V_bias)';
            end

            if show_title == true
                title(set_title,'FontSize', 18);
            end

            if show_max == true
                hold on;
                [max_rho,row,col] = get_max_rho(obj);
                x_max = x(col);
                y_max = y(row);
                plot(x_max, y_max, 'x', 'Color', 'k', 'MarkerSize', 18, 'LineWidth', 3);
                text(x_max+0.2*x(1), y_max-0.2*y(1), sprintf(['max rho = %.2f uC/MV^2 \n' ...
                    'Ec = %.2f MV/cm \n' ...
                    'E_{bias} = %.2f MV/cm'], max_rho,x_max,y_max), ...
                    'Color', 'w', 'FontSize', 18,Interpreter='tex');
                hold off;
            end

            colormap(color);
            % Set equal aspect ratio to make the colormap square
            pbaspect([1 1 1]);
            % Reverse the Y-axis
            set(gca,'YDir','normal','FontSize', 12);
            % Ensure all integers within xlim and ylim are shown
            xticks(floor(min(xlim)):ceil(max(xlim)));
            yticks(floor(min(ylim)):ceil(max(ylim)));
        end

        function [max_rho,row,col] = get_max_rho(obj)

            rho = abs(obj.Rho_E_trans);
            [max_rho, linear_idx] = max(rho(:));
            [row, col] = ind2sub(size(rho), linear_idx);

        end

        function [Ec,Vc] = get_Ec_Vc(obj)
            rho = abs(obj.Rho_E_trans);
            thickness = obj.Thickness;
            vc = obj.Vc_grid;
            [~, linear_idx] = max(rho(:));
            [row, col] = ind2sub(size(rho), linear_idx);
            Vc = vc(row, col); % unit V
            Ec = Vc/thickness * 1e1; % unit MV/cm
        end

        function [Eb,Vb] = get_Eb_Vb(obj)
            rho = abs(obj.Rho_E_trans);
            thickness = obj.Thickness;
            vb = obj.Vb_grid;
            [~, linear_idx] = max(rho(:));
            [row, col] = ind2sub(size(rho), linear_idx);
            Vb = vb(row, col); % unit V
            Eb = Vb/thickness * 1e1; % unit MV/cm
        end

        function [Pr_p,Pr_n] = get_2Pr(obj,current_compensation)
            if ~exist("current_compensation","var")
                current_compensation = true;
            end
            % 提取最后一个VIT loop
            V_end = obj.Voltage_split{end};
            I_end = obj.Current_split{end};
            T_end = obj.Time_split{end}; % unit s
            T_end = T_end - T_end(1);
            area = obj.Area;
            I_end = I_end * 1e6; % unit uA
            Q = cumtrapz(T_end, I_end); % unit uC

            if current_compensation == true
                I_bias_step = 1e-6 * max(abs(I_end));

                bias_val = Q(1) - Q(end);
                n_bias = 1;
                condition_met = abs(Q(1) - Q(end)) / max(abs(Q)) < 0.01; % 闭口偏差小于1%
                while ~condition_met
                    if bias_val > 0
                        Q_sign = 1;
                    else
                        Q_sign = -1;
                    end

                    I = I_end + Q_sign * I_bias_step * n_bias;
                    Q = cumtrapz(T_end,I);

                    bias_val = Q(1) - Q(end);
                    n_bias = n_bias + 1;
                    condition_met = abs(Q(1) - Q(end)) / max(abs(Q)) < 0.01; % 闭口偏差小于1%
                end

            end

            P = Q/area; % unit uC/cm^2
            P = P - 0.5 * (max(P)+min(P));
            % figure;
            % plot(V_end,P);
            idx = find(V_end .* [V_end(end);V_end(1:end-1)] < 0);
            Pr = 0.5 * (P(idx) + P(idx+1));
            if Pr(1) > 0
                Pr_p = Pr(1);
                Pr_n = Pr(2);
            else
                Pr_p = Pr(2);
                Pr_n = Pr(1);
            end

        end

    end
end