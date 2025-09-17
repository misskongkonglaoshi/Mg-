classdef ParticleState < handle
    % ParticleState 管理反应颗粒的状态 (纯数据容器)
    
    properties
        % 几何参数
        r_c                     % 金属核心半径 (m)
        r_p                     % 总半径 (包括氧化层) (m)
        oxide_thickness         % 氧化层厚度 (m)
        initial_oxide_thickness % 初始氧化层厚度 (m)
        
        % 热力学参数
        T_p                    % 颗粒温度 (氧化层表面温度) (K)
        T_c                    % 金属核心温度 (K)
        energy_accumulated     % 熔融过程中累积的能量 (J)
        melting_energy_needed  % 完全熔融所需的能量 (J)
        
        % 反应参数
        reaction_area_factor   % 有效反应面积系数
        is_melting = false    % 是否处于熔融过程
        is_reacting = false   % 是否发生化学反应
        
        % 质量参数
        m_mg               % 金属Mg质量 (kg)
        m_mgo              % MgO质量 (kg)
        m_c                % C质量 (kg)
        melted_fraction       % 熔融金属的质量分数
        
        % 热量参数
        heat_convection = 0    % 对流换热 (W)
        heat_radiation = 0     % 辐射换热 (W)
        heat_reaction = 0      % 反应热 (W)
        heat_reaction_surface = 0      % 反应热 (W)
        heat_total = 0         % 总热量 (W)
    end
    
    properties (Access = private)
        % 参数对象
        params
    end
    
    methods
        function obj = ParticleState(params)
            % 构造函数：初始化颗粒状态
            %
            % 输入:
            %   params: 参数结构体
            
            % 从params初始化几何参数
            if nargin > 0
                % 保存参数
                obj.params = params;
                
                % 初始半径
                initial_radius = params.initial_diameter / 2;
                %fprintf('\n 初始直径 %6e \n',params.initial_diameter);
                % 设置金属核心半径（考虑初始氧化层）
                %if isfield(params, 'initial_oxide_thickness')
                %    fprintf('\n 初始存在氧化膜:\n');
                    obj.initial_oxide_thickness = params.initial_oxide_thickness;
                    obj.r_c = initial_radius - obj.initial_oxide_thickness;
                %else
                %    fprintf('\n 初始wu氧化膜:\n');
                %    obj.initial_oxide_thickness = 0;
                %    obj.r_c = initial_radius;
                %end
                
                % 设置总半径
                obj.r_p = initial_radius;
                
                % 设置氧化层厚度
                obj.oxide_thickness = obj.initial_oxide_thickness;
                
                % 设置温度
                obj.T_p = params.initial_temperature;
                obj.T_c = obj.T_p; % 初始化金属核心温度
                
                % 设置反应面积因子
                if isfield(params, 'reaction_area_factor')
                    obj.reaction_area_factor = params.reaction_area_factor;
                else
                    obj.reaction_area_factor = 1.0;
                end
                
                % 计算初始Mg质量（基于金属核心体积）
                v_mg = (4/3) * pi * obj.r_c^3;
                obj.m_mg = v_mg * obj.params.materials.Mg.density_low;
                
                % 计算初始MgO质量（如果有初始氧化层）
                if obj.initial_oxide_thickness > 0
                    v_total = (4/3) * pi * obj.r_p^3;
                    v_oxide = v_total - v_mg;
                    % 假设初始氧化层全部是MgO
                    obj.m_mgo = v_oxide * obj.params.materials.MgO.density;
                else
                    obj.m_mgo = 0;
                end
                %fprintf('\n 初始mgo质量 %6e \n',obj.initial_oxide_thickness);
                % 初始化其他参数
                obj.m_c = 0;  % 初始无碳
                obj.melted_fraction = 0;
                obj.energy_accumulated = 0;
                obj.is_melting = false;
                obj.is_reacting = false;
            else
                % 设置默认值
                obj.r_c = 0.0;
                obj.r_p = 0.0;
                obj.oxide_thickness = 0.0;
                obj.initial_oxide_thickness = 0.0;
                obj.T_p = 300.0;
                obj.T_c = 300.0; % 默认金属核心温度
                obj.reaction_area_factor = 1.0;
                obj.m_mg = 0.0;
                obj.m_mgo = 0.0;
                obj.m_c = 0.0;
                obj.melted_fraction = 0.0;
                obj.energy_accumulated = 0;
                obj.is_melting = false;
                obj.is_reacting = false;
            end
        end
        
        function update_reaction_area_factor(obj)
            % 根据氧化层厚度更新有效反应面积系数
            if obj.oxide_thickness <= obj.initial_oxide_thickness
                obj.reaction_area_factor = 1.0;
            else
                obj.reaction_area_factor = obj.initial_oxide_thickness / obj.oxide_thickness;
            end
        end
        
        function area = calculate_effective_area(obj)
            % 计算有效反应面积
            % Returns:
            %   area: 有效反应面积 (m²)
            
            % 计算总表面积
            surface_area = 4 * pi * obj.r_p^2;
            
            % 应用有效反应面积系数
            area = surface_area * obj.reaction_area_factor;
        end
        
        function start_melting(obj)
            % 开始熔融过程
            obj.is_melting = true;
            obj.energy_accumulated = 0;
            obj.melting_energy_needed = obj.m_mg * obj.params.materials.Mg.latent_heat;
        end
        
        function [delta_T, is_melting_complete] = process_heat(obj, heat_amount, physicalModel)
            % 处理热量，返回温度变化
            % 输入:
            %   heat_amount: 热量 (J)
            %   physicalModel: PhysicalModel 实例
            % 输出:
            %   delta_T: 温度变化 (K)
            %   is_melting_complete: 是否完成熔融
            
            if obj.is_melting
                % 处于熔融过程中
                obj.energy_accumulated = obj.energy_accumulated + heat_amount;
                
                if obj.energy_accumulated >= obj.melting_energy_needed
                    % 完成熔融
                    obj.is_melting = false;
                    obj.melted_fraction = 1.0;
                    excess_energy = obj.energy_accumulated - obj.melting_energy_needed;
                    obj.energy_accumulated = 0;
                    
                    % 计算剩余能量导致的温度变化
                    total_heat_capacity = physicalModel.calculate_total_heat_capacity(obj);
                    delta_T = excess_energy / total_heat_capacity;
                    is_melting_complete = true;
                else
                    % 继续熔融
                    obj.melted_fraction = obj.energy_accumulated / obj.melting_energy_needed;
                    delta_T = 0;
                    is_melting_complete = false;
                end
            else
                % 正常加热/冷却
                total_heat_capacity = physicalModel.calculate_total_heat_capacity(obj);
                delta_T = heat_amount / total_heat_capacity;
                is_melting_complete = false;
            end
        end
        
        function update_geometry_from_mass(obj)
            % 根据当前质量更新几何参数（半径、厚度）
            
            % 1. 计算Mg核心的体积和半径
            v_mg = obj.m_mg / obj.params.materials.Mg.density;
            if v_mg < 1e-20 % 防止负数开根
                v_mg = 0;
            end
            obj.r_c = (3 * v_mg / (4 * pi))^(1/3);

            % 2. 计算产物层 (MgO+C) 的体积
            v_mgo = obj.m_mgo / obj.params.materials.MgO.density;
            v_c = obj.m_c / obj.params.materials.C.density;
            v_product_shell = v_mgo + v_c;
            
            % 3. 计算总的颗粒体积（核心 + 产物壳层）
            v_total = v_mg + v_product_shell;
            
            % 4. 计算总的颗粒半径
            obj.r_p = (3 * v_total / (4 * pi))^(1/3);

            % 5. 计算产物层厚度
            obj.oxide_thickness = obj.r_p - obj.r_c; % 确保不为负值
            
            % 调试输出
            fprintf('DEBUG: m_mg=%.2e, m_mgo=%.2e, r_c=%.2e, r_p=%.2e, oxide_thickness=%.2e\n', ...
                    obj.m_mg, obj.m_mgo, obj.r_c, obj.r_p, obj.oxide_thickness);
        end
        
        function new_obj = copy(obj)
            % 创建当前对象的一个深拷贝
            new_obj = ParticleState(obj.params); % 调用构造函数以确保所有属性被初始化
            
            % 复制所有公共属性的值
            props = properties(obj);
            for i = 1:length(props)
                prop_name = props{i};
                new_obj.(prop_name) = obj.(prop_name);
            end
            
            % 确保热量字段被复制
            if isfield(obj, 'heat_convection')
                new_obj.heat_convection = obj.heat_convection;
                new_obj.heat_radiation = obj.heat_radiation;
                new_obj.heat_reaction = obj.heat_reaction;
                new_obj.heat_reaction_surface = obj.heat_reaction_surface;
                new_obj.heat_total = obj.heat_total;
            else
                new_obj.heat_convection = 0;
                new_obj.heat_radiation = 0;
                new_obj.heat_reaction = 0;
                new_obj.heat_reaction_surface = 0;
                new_obj.heat_total = 0;
            end
        end
    end
end 