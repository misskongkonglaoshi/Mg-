class SimulationParameters:
    def __init__(self):
        # 颗粒参数
        self.particle_diameter = 2.0e-4  # 初始颗粒直径 (m)
        self.initial_temperature = 300.0  # 初始温度 (K)
        self.initial_oxide_ratio = 0.02   # 初始氧化层厚度占颗粒半径的比例
        
        # 物质密度 (kg/m³)
        self.rho_mg = 1740    # Mg
        self.rho_mgo = 3580   # MgO
        self.rho_c = 2260     # C (石墨)
        
        # 环境参数
        self.ambient_temperature = 1200.0  # 环境温度 (K)
        self.ambient_pressure = 1.0e5      # 环境压力 (Pa)
        self.gas_velocity = 0.1            # 气体流速 (m/s)
        
        # 气体组成 (摩尔分数)
        self.ambient_gas_composition = {
            "CO2": 1.0,  # 纯CO2环境
            "O2": 0.0,
            "N2": 0.0
        }
        
        # 反应动力学参数
        self.A_pre = 1.376e-4  # 指前因子 (s^-1)
        self.E_a = 132.4e3     # 活化能 (J/mol)
        self.T_melting = 923   # Mg熔点 (K)
        self.T_ignition = 1366 # Mg点火温度 (K)
        
        # 时间步长设置
        self.time_step = 1.0e-6  # 时间步长 (s)
        self.total_time = 0.1    # 总模拟时间 (s)
        
        # 输出控制
        self.output_interval = 100  # 每隔多少步输出一次结果
        self.save_results = True    # 是否保存结果
        self.plot_results = True    # 是否绘制结果
        
        # 结果记录参数
        self.record_variables = [
            "time",              # 时间 (s)
            "temperature",       # 温度 (K)
            "metal_radius",      # 金属核心半径 (m)
            "total_radius",      # 总半径 (m)
            "oxide_thickness",   # 氧化层厚度 (m)
            "mass_mg",          # 金属Mg质量 (kg)
            "mass_mgo",         # MgO质量 (kg)
            "mass_c",           # C质量 (kg)
            "reaction_rate",     # 反应速率 (kg/s)
            "q_conv",           # 对流换热 (W)
            "q_rad",            # 辐射换热 (W)
            "q_reaction"        # 反应热 (W)
        ]