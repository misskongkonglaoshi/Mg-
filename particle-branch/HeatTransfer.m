classdef HeatTransfer < handle
    % HeatTransfer 处理颗粒的热传递计算
    
    properties
        params  % 参数对象
    end
    
    methods
        function obj = HeatTransfer(params)
            % 构造函数：初始化热传递模型
            %
            % 输入:
            %   params: 参数对象
            
            obj.params = params;
        end
        
        function q_conv = calculateConvection(obj, T_particle, T_ambient, diameter)
            % 计算对流换热
            %
            % 输入:
            %   T_particle: 颗粒温度 (K)
            %   T_ambient: 环境温度 (K)
            %   diameter: 颗粒直径 (m)
            %
            % 输出:
            %   q_conv: 对流换热热量 (W)
            
            % 计算表面积
            surface_area = pi * diameter^2;
            
            % 计算对流换热 (W)
            q_conv = obj.params.h_conv * surface_area * (T_ambient - T_particle);
        end
        
        function q_rad = calculateRadiation(obj, T_particle, T_ambient, diameter)
            % 计算辐射换热
            %
            % 输入:
            %   T_particle: 颗粒温度 (K)
            %   T_ambient: 环境温度 (K)
            %   diameter: 颗粒直径 (m)
            %
            % 输出:
            %   q_rad: 辐射换热热量 (W)
            
            % 计算表面积
            surface_area = pi * diameter^2;
            
            % 计算辐射换热 (W)
            q_rad = obj.params.emissivity * obj.params.sigma * surface_area * ...
                (T_ambient^4 - T_particle^4);
        end
        
        function q_total = calculateTotalHeatTransfer(obj, T_particle, T_ambient, diameter, oxide_thickness)
            % 计算总热传递，考虑氧化层的影响
            %
            % 输入:
            %   T_particle: 颗粒温度 (K)
            %   T_ambient: 环境温度 (K)
            %   diameter: 颗粒直径 (m)
            %   oxide_thickness: 氧化层厚度 (m)，可选参数
            %
            % 输出:
            %   q_total: 总热传递 (W)
            
            % 计算对流和辐射热传递
            q_conv = obj.calculateConvection(T_particle, T_ambient, diameter);
            q_rad = obj.calculateRadiation(T_particle, T_ambient, diameter);
            
            % 基本热传递
            q_total = q_conv + q_rad;
            
            % 如果提供了氧化层厚度，考虑其影响
            if nargin > 4 && oxide_thickness > 0
                % 考虑氧化层的热阻
                % 热阻 = 厚度 / (导热系数 * 面积)
                surface_area = pi * diameter^2;
                R_oxide = oxide_thickness / (obj.params.materials.MgO.thermal_conductivity * surface_area);
                
                % 计算通过氧化层的最大热流
                delta_T = abs(T_ambient - T_particle);
                q_max = delta_T / R_oxide;
                
                % 取较小值作为实际热传递
                q_total = min(q_total, q_max);
            end
        end
    end
end 