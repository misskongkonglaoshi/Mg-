classdef ThermoReader < handle
    % ThermoReader 读取CHEMKIN-III格式的热力学数据文件并计算性质
    
    properties
        thermo_file
        params
        species_data
        Ru = 8.31446261815324; % J/(mol·K)
        warned_species = {}; % V31: 新增，用于记录已警告过的缺失物种
    end
    
    methods
        function obj = ThermoReader(thermo_file, params)
            if nargin > 0
                obj.thermo_file = thermo_file;
                if nargin > 1
                    obj.params = params;
                end
                obj.species_data = struct();
                obj.read_thermo_file();
            end
        end
        
        function read_thermo_file(obj)
            fid = fopen(obj.thermo_file, 'r');
            if fid == -1
                error('无法打开文件: %s', obj.thermo_file);
            end
            
            % 跳过文件头直到THERMO关键字
            line = fgetl(fid);
            while ischar(line) && ~strncmpi(strtrim(line), 'THERMO', 6)
                line = fgetl(fid);
            end
            
            % 读取并忽略全局温度行
            fgetl(fid);
            
            % 循环读取物种数据
            while true
                line1 = fgetl(fid);
                if ~ischar(line1) || strncmpi(strtrim(line1), 'END', 3)
                    break;
                end
                
                % 跳过文件中的空行或纯空格行
                if isempty(strtrim(line1))
                    continue;
                end
                
                line2 = fgetl(fid);
                line3 = fgetl(fid);
                line4 = fgetl(fid);

                try
                    % --- 开始修改: 使用sprintf填充所有系数行以确保长度 ---
                    line2 = sprintf('%-75s', line2);
                    line3 = sprintf('%-75s', line3);
                    line4 = sprintf('%-75s', line4);
                    % --- 结束修改 ---

                    % 解析第一行
                    species_name = strtrim(line1(1:18));
                    % --- 开始修改: 扩大温度读取范围 ---
                    temp_ranges_str = sscanf(line1(46:73), '%f');
                    % --- 结束修改 ---
                    T_low = temp_ranges_str(1);
                    T_high = temp_ranges_str(2);
                    T_mid = temp_ranges_str(3);

                    % --- 开始修改: 使用健壮的固定宽度解析方法 ---
                    coeffs_str_cell = {};
                    % 第2行: 5个系数
                    for i=1:5
                        coeffs_str_cell{end+1} = line2((i-1)*15+1 : i*15);
                    end
                    % 第3行: 5个系数
                    for i=1:5
                        coeffs_str_cell{end+1} = line3((i-1)*15+1 : i*15);
                    end
                    % 第4行: 4个系数
                    for i=1:4
                        coeffs_str_cell{end+1} = line4((i-1)*15+1 : i*15);
                    end
                    coeffs_vals = str2double(coeffs_str_cell);
                    % --- 结束修改 ---
                    
                    if length(coeffs_vals) ~= 14
                        warning('物种 %s 的系数数量不为14，跳过。', species_name);
                        continue;
                    end
                    
                    % NASA 7系数多项式
                    high_T_coeffs = coeffs_vals(1:7);
                    low_T_coeffs = coeffs_vals(8:14);
                    
                    obj.species_data.(species_name) = struct(...
                        'T_low', T_low, 'T_high', T_high, 'T_mid', T_mid, ...
                        'high_T_coeffs', high_T_coeffs, ...
                        'low_T_coeffs', low_T_coeffs);
                catch ME
                    warning('解析物种 %s 时出错: %s', species_name, ME.message);
                end
            end
            
            fclose(fid);
        end
        
        function coeffs = get_coeffs(obj, species, T)
            if ~isfield(obj.species_data, species)
                error('缺少物种 %s 的热力学数据', species);
            end
            data = obj.species_data.(species);
            if T > data.T_mid
                coeffs = data.high_T_coeffs;
            else
                coeffs = data.low_T_coeffs;
            end
        end

        function Cp = calculate_Cp(obj, species, T)
            % 计算定压比热容 (J/mol·K)
            a = obj.get_coeffs(species, T);
            T_vec = [1; T; T.^2; T.^3; T.^4];
            Cp = (a(1:5) * T_vec) * obj.Ru;
        end
        
        function H = calculate_H(obj, species, T)
            % 计算焓 (J/mol)
            a = obj.get_coeffs(species, T);
            T_vec = [1; T/2; T.^2/3; T.^3/4; T.^4/5; 1./T];
            H = (a(1:6) * T_vec) * obj.Ru * T;
        end
        
        function S = calculate_S(obj, species, T)
            % 计算熵 (J/mol·K)
            a = obj.get_coeffs(species, T);
            T_vec = [log(T); T; T.^2/2; T.^3/3; T.^4/4];
            S = (a(1:5) * T_vec + a(7)) * obj.Ru;
        end
        
        function p_sat = get_p_sat(obj, T_kelvin)
            if isempty(obj.params)
                 error('需要params对象来计算饱和蒸气压');
            end
            T_celsius = T_kelvin - 273.15;
            A = obj.params.antoine_coeffs_mg.A;
            B = obj.params.antoine_coeffs_mg.B;
            C = obj.params.antoine_coeffs_mg.C;
            log10_p_mmHg = A - B / (T_celsius + C);
            p_mmHg = 10^log10_p_mmHg;
            p_sat = p_mmHg * 133.322;
        end
        
        function has_data = has_species(obj, species)
            has_data = isfield(obj.species_data, species);
        end

        function dH = get_reaction_enthalpy(obj, reaction_string, T)
            % 计算化学反应的焓变 (J/mol of reaction)
            %
            % 输入:
            %   reaction_string: 形如 "A+2B=C+D" 的字符串
            %   T: 温度 (K)
            % 输出:
            %   dH: 反应焓变 (J/mol)

            % 1. 解析反应物和产物
            parts = strsplit(reaction_string, '=');
            if numel(parts) ~= 2
                error('ThermoReader:InvalidReactionFormat', ...
                      '反应式 "%s" 格式无效，必须包含一个 "="。', reaction_string);
            end
            
            reactant_str = parts{1};
            product_str = parts{2};

            % 2. 计算产物总焓
            H_products = obj.calculate_side_enthalpy(product_str, T);

            % 3. 计算反应物总焓
            H_reactants = obj.calculate_side_enthalpy(reactant_str, T);
            
            % 4. 返回反应焓变
            dH = H_products - H_reactants;
        end
        
        function H_total = calculate_side_enthalpy(obj, side_string, T)
            % Helper: 计算反应式一侧的总焓
            terms = strsplit(side_string, '+');
            H_total = 0;
            for i = 1:numel(terms)
                term = strtrim(terms{i});
                [coeff, species_name] = obj.parse_term(term);
                
                % 将用户友好的名称映射到thermo.dat中的键
                [thermo_key, is_solid] = obj.map_species_name(species_name);

                % 检查物种是否存在
                if obj.has_species(thermo_key)
                    % 优先使用thermo.dat中的详细数据
                    H_total = H_total + coeff * obj.calculate_H(thermo_key, T);
                elseif is_solid && isfield(obj.params.materials, thermo_key)
                    % 如果是固相且在thermo.dat中找不到, 则回退到params.m中的Hf298
                    % 这是一个近似, 忽略了固相焓随温度的变化
                    
                    % V31: 增加只警告一次的逻辑
                    if ~ismember(species_name, obj.warned_species)
                        fprintf('注意: 在thermo.dat中未找到固相 "%s", 已回退使用params.m中定义的Hf298。\n', species_name);
                        obj.warned_species{end+1} = species_name;
                    end

                    H_total = H_total + coeff * obj.params.materials.(thermo_key).Hf298;
                else
                    error('ThermoReader:SpeciesNotFound', ...
                          '在解析 "%s" 时，无法在任何数据源中找到物种 "%s" (映射为 "%s")', ...
                          term, species_name, thermo_key);
                end
            end
        end

        function [coeff, species_name] = parse_term(~, term)
            % Helper: 解析 "2H2O" 这样的项
            % 使用正则表达式查找前导数字 (系数) 和物种名称
            tokens = regexp(term, '^\s*(\d*\.?\d*)\s*([\w\(\)]+)\s*$', 'tokens');
            
            if isempty(tokens)
                error('ThermoReader:InvalidTermFormat', '无法解析反应项: "%s"', term);
            end
            
            coeff_str = tokens{1}{1};
            species_name = tokens{1}{2};
            
            if isempty(coeff_str)
                coeff = 1;
            else
                coeff = str2double(coeff_str);
            end
        end

        function [thermo_key, is_solid] = map_species_name(~, user_name)
            % Helper: 将用户输入的名称映射到 thermo.dat 中的键, 并判断是否为固相
            % 'Mg(g)' -> 'Mg', 'MgO(s)' -> 'MgO', 'C(s)' -> 'C', etc.
            
            is_solid = contains(user_name, '(s)', 'IgnoreCase', true);
            
            % 使用正则表达式移除末尾的相态标识, 例如 (s), (g), (l)
            % 这会保留物种名称本身的大小写, 例如 MgO
            key = regexprep(user_name, '\s*\([sglSGL]\)\s*$', '');
            
            thermo_key = strtrim(key);
        end

        function value = get_property(obj, species, property_name)
            % 通用属性获取函数
            % 优先从thermo.dat计算，然后回退到params.m
            
            % 尝试从CHEMKIN数据计算
            if obj.has_species(species)
                switch lower(property_name)
                    case 'hf298' % 298.15K下的标准生成焓
                        value = obj.calculate_H(species, 298.15);
                        return;
                    case 'cp'
                        % 注意: Cp是温度依赖的, 这里只为示例，实际应调用calculate_Cp
                        % value = obj.calculate_Cp(species, T); 
                        % 此处需要一个温度参数, 因此不直接处理
                end
            end
            
            % 如果无法计算, 尝试从params回退
            if ~isempty(obj.params) && isfield(obj.params.materials, species) ...
               && isfield(obj.params.materials.(species), property_name)
                value = obj.params.materials.(species).(property_name);
                return;
            end
            
            % 如果都找不到，则报错
            error('ThermoReader:PropertyNotFound', ...
                  '无法在thermo.dat或Parameters.m中找到物种 "%s" 的属性 "%s"。', ...
                  species, property_name);
        end
    end
end 