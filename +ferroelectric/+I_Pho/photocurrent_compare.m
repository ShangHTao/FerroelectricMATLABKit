classdef photocurrent_compare
    properties
        File_path = '';
        Sample_name = '';
        Setup_title = {};
        CVS_file_path = {};
        Voltage_Prepolar = [];
        Photocurrent_objs = {};
        N_data = 0;
    end

    methods
        function obj = photocurrent_compare(file_path,sample_name)
            import ferroelectric.I_Pho.*;

            obj.File_path = file_path;
            obj.Sample_name = sample_name;
            cvs_files = dir(fullfile(file_path,'I_V-t Sampling*.csv'));

            N = numel(cvs_files);
            objs = cell(N,1);
            voltage_prepolar = zeros(N,1);

            for i = 1:N
                token = regexp(cvs_files(i).name, 'after\s*([+-]?\d+\.?\d*)V', 'tokens', 'once');
                voltage_prepolar(i) = str2double(token{1});
                cvs_file_path = [cvs_files(i).folder filesep cvs_files(i).name];
                objs{i} = photocurrent_single(cvs_file_path);
                obj.Setup_title{i} = objs{i}.Setup_title;
            end
            obj.N_data = N;
            obj.Photocurrent_objs = objs;
            obj.Voltage_Prepolar = voltage_prepolar;
        end

        function obj = plot_compare(obj)
            N = obj.N_data;
            v_prepo = obj.Voltage_Prepolar;
            objs = obj.Photocurrent_objs;
            samplename = obj.Sample_name; 
            legend_str = cell(1,2*N);
            color = [1 0 0 % r
                0 1 0 % g
                0 0 1];% b

            for i = 1:N
                t_s = objs{i}.Time_raw;
                % v = obj.Voltage_raw;
                c_pA = objs{i}.Current_raw*1e12;
                plot(t_s, c_pA,'LineWidth', 1.5, 'LineStyle', '-', 'Marker', 'none', 'MarkerSize', 8, 'Color', [color(i,:) 0.3]);
                hold on;
                if ~exist('filter_window','var')
                    filter_window = round(0.1*numel(c_pA));
                end
                I_filtered = movmean(c_pA, filter_window);
                plot(t_s, I_filtered,'LineWidth', 2.5, 'LineStyle', '-', 'Marker', 'none', 'MarkerSize', 8, 'Color', color(i,:));
                hold on;
                legend_str{2*i-1} = ['V_Prepo:',num2str(v_prepo(i)), 'V raw data'];
                legend_str{2*i} = ['V_Prepo:',num2str(v_prepo(i)), 'V filter data'];

                % disp([obj.Sample_name,' ' ,num2str(v_prepo(i)),'V mean current = ',  num2str(mean(c_pA)), ' pA']);
            end
            legend(legend_str,Interpreter="none",Location = "best");
            xlabel('Time$(s)$', 'FontSize', 18,'Interpreter','latex');
            ylabel('$I_{pho}(pA)$', 'FontSize', 18,'Interpreter','latex');
            xlim([min(t_s),max(t_s)]);
            ymax = get_absmax_current(obj);
            ylim([-1.3*ymax*1e12,1.3*ymax*1e12]);
            title(['$I_{pho}$ compare ' samplename],'FontSize', 18,'Interpreter','latex');
        end

        function I_absmax = get_absmax_current(obj)
            N = obj.N_data;
            photocurrent_objs = obj.Photocurrent_objs;
            I_absmax_n = zeros(N,1);
            for i = 1:N
                I_absmax_n(i) = photocurrent_objs{i}.Current_absmax;
            end
            I_absmax = max(abs(I_absmax_n));
        end

        function T = get_current_table(obj,filter_window)
            N = obj.N_data;
            v_prepo = obj.Voltage_Prepolar;
            objs = obj.Photocurrent_objs;
            samples = strings(N,1);
            current = cell(N,1);
            current_filtered = cell(N,1);
            current_mean = zeros(N,1);
            for i = 1:N
                samples(i) = v_prepo(i);
                c = objs{i}.Current_raw;
                current{i} = c;
                if ~exist('filter_window','var')
                    filter_window = round(0.1*numel(c));
                end
                current_filtered{i} = movmean(c, filter_window);
                current_mean(i) = mean(c);

            end
            T = table(samples,current,current_filtered,current_mean, ...
                'VariableNames',{'Prepolar_v','Current','Current_filtered','Current_mean'});
        end


    end
end