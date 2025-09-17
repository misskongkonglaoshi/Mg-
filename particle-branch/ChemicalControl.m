classdef ChemicalControl < handle
    % ChemicalControl 处理化学反应动力学控制
    
    properties
        % 物理常数
        Ru = 8314.5;  % 普适气体常数 (J/(kmol·K))
        
        % 反应相关参数
        h_reaction = -601.6e3;  % Mg + CO2 -> MgO + C 的反应焓 (J/mol)
        M_mg = 24.305;  % Mg的摩尔质量 (g/mol)
        MgO = 40.304;
        M_co2 = 44.01;  % CO2的摩尔质量 (g/mol)
        A_pre = 1.376e-4;  % 指前因子 (s-1)
        E_a = 132.4e3;  % 活化能 (J/mol)
        
        % 反应计量系数
        stoichiometry
        
        % 热力学数据读取器
        thermo_reader
        
        % 颗粒状态
        particle_state
    end
    
    methods
        function obj = ChemicalControl(particle_state, params)
            % 构造函数：初始化化学反应控制模型
            %
            % 输入:
            %   particle_state: ParticleState 实例
            %   params: 其他参数
            
            % 设置颗粒状态
            obj.particle_state = particle_state;
            
            % 初始化反应计量系数
            obj.stoichiometry = struct();
            obj.stoichiometry.reactants = struct('Mg', 2, 'CO2', 1);
            obj.stoichiometry.products = struct('MgO', 2, 'C', 1);
            
            % 初始化热力学数据读取器
            thermo_file = fullfile(fileparts(mfilename('fullpath')), 'thermo.dat');
            obj.thermo_reader = ThermoReader(thermo_file);
        end
        
        function coeffs = get_coeffs(obj, species, T)
            % 获取对应温度下的NASA系数
            coeffs = obj.thermo_reader.get_coeffs(species, T);
        end
        
        function G = calculate_gibbs(obj, species, T)
            % 计算特定温度下的摩尔吉布斯自由能 (J/mol)
            G = obj.calculate_enthalpy(species, T) - T * obj.calculate_entropy(species, T);
        end
        
        function dH = calculate_reaction_enthalpy(obj, T)
            % 计算反应焓变，考虑化学计量数
            % 计算产物总焓
            H_products = 0;
            product_species = fieldnames(obj.stoichiometry.products);
            for i = 1:length(product_species)
                species = product_species{i};
                coeff = obj.stoichiometry.products.(species);
                H_products = H_products + coeff * obj.calculate_enthalpy(species, T);
            end
            
            % 计算反应物总焓
            H_reactants = 0;
            reactant_species = fieldnames(obj.stoichiometry.reactants);
            for i = 1:length(reactant_species)
                species = reactant_species{i};
                coeff = obj.stoichiometry.reactants.(species);
                H_reactants = H_reactants + coeff * obj.calculate_enthalpy(species, T);
            end
            
            % 返回反应焓变
            dH = H_products - H_reactants;
        end
        
        function dG = calculate_reaction_gibbs(obj, T)
            % 计算反应吉布斯自由能变，考虑化学计量数
            % 计算产物总吉布斯自由能
            G_products = 0;
            product_species = fieldnames(obj.stoichiometry.products);
            for i = 1:length(product_species)
                species = product_species{i};
                coeff = obj.stoichiometry.products.(species);
                G_products = G_products + coeff * obj.calculate_gibbs(species, T);
            end
            
            % 计算反应物总吉布斯自由能
            G_reactants = 0;
            reactant_species = fieldnames(obj.stoichiometry.reactants);
            for i = 1:length(reactant_species)
                species = reactant_species{i};
                coeff = obj.stoichiometry.reactants.(species);
                G_reactants = G_reactants + coeff * obj.calculate_gibbs(species, T);
            end
            
            % 返回反应吉布斯自由能变
            dG = G_products - G_reactants;
        end
        
        function rate = calculate_reaction_rate(obj, T, p_co2, oxide_thickness)
            % 计算化学反应速率，考虑氧化层影响
            %
            % 输入:
            %   T: 温度 (K)
            %   p_co2: CO2分压 (Pa)
            %   oxide_thickness: 氧化层厚度 (m)，可选参数
            %
            % 输出:
            %   rate: 反应速率 (kg/s)
            
            % 计算基本反应速率常数
            k = obj.calculate_rate_constant(T);
            
            % 如果颗粒状态已设置，获取有效反应面积；否则使用默认值
            if ~isempty(obj.particle_state)
                effective_area = obj.particle_state.calculate_effective_area();
            else
                % 假设默认有效面积
                effective_area = 1e-8;  % 默认值
            end
            
            % 计算基本反应速率
            rate = k * p_co2 * effective_area;
            
            % 如果提供了氧化层厚度，考虑其影响
            if nargin > 3 && oxide_thickness > 0
                % 氧化层阻碍效应（简化模型：指数衰减）
                oxide_factor = exp(-oxide_thickness / 1e-6);  % 使用1µm作为特征长度
                rate = rate * oxide_factor;
            end
        end
        
        function k = calculate_rate_constant(obj, T)
            % 将气体常数从 J/(kmol·K) 转换为 J/(mol·K)
            R = obj.Ru / 1000;
            
            % 计算反应速率常数
            k = obj.A_pre * exp(-obj.E_a / (R * T));
        end
    end
end 