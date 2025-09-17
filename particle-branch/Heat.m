classdef Heat < handle
    % Heat 处理颗粒的热传递计算
    
    properties
        % 物理常数
        kb = 1.380649e-23;  % 玻尔兹曼常数 (J/K)
        sigma = 5.67e-8;    % Stefan-Boltzmann常数 (W/(m²·K⁴))
        
        % CO2气体参数
        d_co2 = 3.3e-10;    % CO2分子直径 (m)
        k_co2 = 0.0146;     % CO2导热系数 (W/(m·K)) at 273K
        cp_co2 = 844;       % CO2比热容 (J/(kg·K))
        mu_co2 = 1.47e-5;   % CO2动力粘度 (Pa·s)
        
        % 辐射参数
        emissivity_params = struct(...
            'a', 0.7,       % 基础发射率
            'b', 0.0001     % 温度系数
        );
        
        % 颗粒状态
        particle_state
        
        % 环境参数
        params
    end
    
    methods
        function obj = Heat(particle_state, params)
            % 构造函数：初始化热传递模型
            %
            % 输入:
            %   particle_state: ParticleState实例
            %   params: 环境参数
            
            obj.particle_state = particle_state;
            obj.params = params;
        end
        
        function mass_fractions = convert_mole_to_mass_fraction(obj, gas_components)
            % 将摩尔分数转换为质量分数
            % 输入:
            %   gas_components: 气体组分结构体
            % 输出:
            %   mass_fractions: 质量分数结构体
            
            % 计算总摩尔质量
            total_mass = 0;
            species = fieldnames(gas_components);
            for i = 1:length(species)
                species_name = species{i};
                total_mass = total_mass + gas_components.(species_name).mole_fraction * ...
                    gas_components.(species_name).M_w;
            end
            
            % 计算质量分数
            mass_fractions = struct();
            for i = 1:length(species)
                species_name = species{i};
                mass_fractions.(species_name) = gas_components.(species_name).mole_fraction * ...
                    gas_components.(species_name).M_w / total_mass;
            end
        end
        
        function mu = calculate_single_viscosity(obj, species, T)
            % 计算单组元气体的粘性系数
            % 输入:
            %   species: 气体组分名称
            %   T: 温度 (K)
            % 输出:
            %   mu: 粘性系数 (Pa·s)
            
            % 获取气体参数
            M_w = obj.params.gas_components.(species).M_w;
            sigma_d = obj.params.gas_components.(species).sigma_d;
            epsilon_min = obj.params.gas_components.(species).epsilon_min;
            
            % 计算约化温度
            T_star = sigma_d * T / epsilon_min;
            
            % 计算碰撞积分（简化模型）
            Omega_u = 1.16145 * T_star^(-0.14874) + 0.52487 * exp(-0.77320 * T_star) + ...
                2.16178 * exp(-2.43787 * T_star);
            
            % 计算粘性系数
            mu = 2.6693e-6 * sqrt(M_w * T) / (sigma_d^2 * Omega_u);
        end
        
        function k = calculate_single_thermal_conductivity(obj, species, T, mu)
            % 计算单组元气体的导热系数
            % 输入:
            %   species: 气体组分名称
            %   T: 温度 (K)
            %   mu: 粘性系数 (Pa·s)
            % 输出:
            %   k: 导热系数 (W/(m·K))
            
            % 获取气体参数
            M_w = obj.params.gas_components.(species).M_w;
            cp = obj.params.gas_components.(species).cp;
            
            % 计算导热系数
            k = 15/4 * obj.params.R_u/M_w * mu * (4/15 * M_w * cp/obj.params.R_u + 1/3);
        end
        
        function [mu_mix, k_mix] = calculate_mixture_properties(obj, T)
            % 计算混合气体的物性参数
            % 输入:
            %   T: 温度 (K)
            % 输出:
            %   mu_mix: 混合气体粘性系数 (Pa·s)
            %   k_mix: 混合气体导热系数 (W/(m·K))
            
            % 获取质量分数
            mass_fractions = obj.convert_mole_to_mass_fraction(obj.params.gas_components);
            
            % 初始化混合物性质
            mu_mix = 0;
            k_mix = 0;
            
            % 计算每种气体的贡献
            species = fieldnames(obj.params.gas_components);
            for i = 1:length(species)
                species_name = species{i};
                
                % 计算单组元气体的物性
                mu = obj.calculate_single_viscosity(species_name, T);
                k = obj.calculate_single_thermal_conductivity(species_name, T, mu);
                
                % 按质量分数加权
                mu_mix = mu_mix + mass_fractions.(species_name) * mu;
                k_mix = k_mix + mass_fractions.(species_name) * k;
            end
        end
        
        function Kn = calculate_knudsen(obj, T)
            % 计算Knudsen数
            % 输入:
            %   T: 温度 (K)
            % 输出:
            %   Kn: Knudsen数
            
            % 计算平均分子直径（使用混合气体的平均值）
            d_avg = 0;
            total_mole_fraction = 0;
            species = fieldnames(obj.params.gas_components);
            for i = 1:length(species)
                species_name = species{i};
                d_avg = d_avg + obj.params.gas_components.(species_name).mole_fraction * ...
                    obj.params.gas_components.(species_name).sigma_d;
                total_mole_fraction = total_mole_fraction + ...
                    obj.params.gas_components.(species_name).mole_fraction;
            end
            d_avg = d_avg / total_mole_fraction;
            
            % 计算平均自由程
            lambda = obj.params.kb * T / (sqrt(2) * pi * d_avg^2 * obj.params.ambient_pressure);
            
            % 计算特征长度（颗粒直径）
            L = 2 * obj.particle_state.total_radius;
            
            % 计算Knudsen数
            Kn = lambda / L;
        end
        
        function h = calculate_heat_transfer_coefficient(obj, T_particle, T_ambient)
            % 根据Knudsen数计算对流换热系数
            % 输入:
            %   T_particle: 颗粒温度 (K)
            %   T_ambient: 环境温度 (K)
            % 输出:
            %   h: 对流换热系数 (W/(m²·K))
            
            % 计算Knudsen数
            Kn = obj.calculate_knudsen(T_particle);
            
            % 计算特征长度
            d = 2 * obj.particle_state.total_radius;
            
            % 计算混合气体的物性参数
            T_film = (T_particle + T_ambient) / 2;
            [mu_mix, k_mix] = obj.calculate_mixture_properties(T_film);
            
            if Kn < 0.01
                % 连续区域：使用传统Nu数计算
                Re = obj.params.ambient_velocity * d * obj.params.ambient_density / mu_mix;
                Pr = mu_mix * obj.params.cp_mix / k_mix;
                Nu = 2 + 0.6 * Re^0.5 * Pr^0.33;
                h = Nu * k_mix / d;
                
            elseif Kn < 0.1
                % 滑移区域：考虑温度跳跃
                beta_t = 2.2;  % 温度跳跃系数
                h = k_mix / (d + 2 * beta_t * obj.calculate_knudsen(T_film) * d);
                
            elseif Kn < 10
                % 过渡区域：使用过渡区域修正
                alpha = 1;  % 热适应系数
                gamma = 1.4;  % 比热比
                Pr = mu_mix * obj.params.cp_mix / k_mix;
                h = alpha * (gamma + 1)/(8 * (gamma - 1)) * ...
                    sqrt(obj.params.R_u * T_film / (2 * pi * obj.params.M_w_mix)) * ...
                    obj.params.ambient_density * obj.params.cp_mix;
                
            else
                % 自由分子区域：使用自由分子理论
                alpha = 1;  % 热适应系数
                h = alpha * obj.params.ambient_density * obj.params.cp_mix * ...
                    sqrt(obj.params.R_u * T_film / (2 * pi * obj.params.M_w_mix)) / 4;
            end
        end
        
        function emissivity = calculate_emissivity(obj, T)
            % 计算温度相关的发射率
            % 输入:
            %   T: 温度 (K)
            % 输出:
            %   emissivity: 发射率
            
            % 使用参数中的固定发射率
            emissivity = obj.params.emissivity;
        end
        
        function q_radiation = calculate_radiation_heat(obj, T_ambient)
            % 计算辐射换热热量
            % 输入:
            %   T_ambient: 环境温度 (K)
            % 输出:
            %   q_radiation: 辐射换热热量 (W)
            
            % 计算表面积
            surface_area = 4 * pi * obj.particle_state.total_radius^2;
            
            % 计算发射率
            emissivity = obj.calculate_emissivity(obj.particle_state.temperature);
            
            % 计算辐射换热 (W)
            q_radiation = obj.params.sigma * emissivity * surface_area * ...
                (T_ambient^4 - obj.particle_state.temperature^4);
        end
        
        function q_convection = calculate_convection_heat(obj, T_ambient)
            % 计算对流换热热量
            % 输入:
            %   T_ambient: 环境温度 (K)
            % 输出:
            %   q_convection: 对流换热热量 (W)
            
            % 计算表面积
            surface_area = 4 * pi * obj.particle_state.total_radius^2;
            
            % 计算对流换热系数
            h = obj.calculate_heat_transfer_coefficient(...
                obj.particle_state.temperature, T_ambient);
            
            % 计算对流换热 (W)
            q_convection = h * surface_area * ...
                (T_ambient - obj.particle_state.temperature);
        end
    end
end
