import numpy as np

class DiffusionControl:
    def __init__(self):
        # 物理常数
        self.Ru = 8314.5  # 普适气体常数 (J/(kmol·K))
        self.sc = 0.6  # Schmidt数
        
        # 气体特性
        self.gas_properties = {
            "O2": {"M_w": 32.00, "sigma_d": 3.467e-10, "epsilon_min": 106.7},
            "CO2": {"M_w": 44.01, "sigma_d": 3.941e-10, "epsilon_min": 195.2},
            "H2O": {"M_w": 18.01528, "sigma_d": 2.641e-10, "epsilon_min": 809.1},
            "N2": {"M_w": 28.0134, "sigma_d": 3.798E-10, "epsilon_min": 71.4},
            "MgO": {"M_w": 40.304, "sigma_d": 3.941e-10, "epsilon_min": 195.2},
            "C": {"M_w": 12.011, "sigma_d": 3.941e-10, "epsilon_min": 195.2}
        }

    def calculate_collision_integral(self, T_star):
        """计算碰撞积分"""
        collision_integral = 1.147 * T_star**(-0.145) + (T_star + 0.5)**(-2)
        return collision_integral

    def calculate_diffusion_coefficient(self, temp, pressure, species):
        """计算扩散系数"""
        if species not in self.gas_properties:
            raise ValueError(f"Species {species} not found in gas properties")
            
        # 获取气体特性
        M_w = self.gas_properties[species]["M_w"]
        sigma_d = self.gas_properties[species]["sigma_d"]
        epsilon_min = self.gas_properties[species]["epsilon_min"]
        
        # 计算扩散系数
        T_star = (sigma_d * temp) / epsilon_min
        Omega_D = self.calculate_collision_integral(T_star)
        
        # 使用Chapman-Enskog理论计算二元扩散系数
        D_ab = 1.8583e-7 * np.sqrt(temp**3 * (1/M_w + 1/28.97)) / (pressure * sigma_d**2 * Omega_D)
        
        return D_ab

    def calculate_sherwood_number(self, particle_diameter, gas_velocity, nu, sc):
        """计算Sherwood数，使用Frossling关联式"""
        # 计算雷诺数
        re = (particle_diameter * gas_velocity) / nu
        # Frossling关联式
        sh = 2 + 0.6 * (re ** 0.5) * (sc ** (1/3))
        return sh

    def compute_diffusion_controlled_rate(self, particle_temp, ambient_temp, particle_diameter,
                                        gas_velocity, ambient_composition, ambient_pressure, 
                                        reaction_area_factor, nu):
        """计算扩散控制下的反应速率"""
        # 计算CO2的扩散系数
        D_co2 = self.calculate_diffusion_coefficient(ambient_temp, ambient_pressure, "CO2")
        
        # 计算Sherwood数
        sh = self.calculate_sherwood_number(particle_diameter, gas_velocity, nu, self.sc)
        
        # 计算质量传递系数
        k_m = sh * D_co2 / particle_diameter
        
        # 计算CO2浓度差
        c_inf = (ambient_pressure * ambient_composition.get("CO2", 0)) / (self.Ru * ambient_temp)
        c_s = 0  # 假设表面CO2浓度为0
        
        # 计算反应面积
        surface_area = np.pi * particle_diameter**2
        effective_area = surface_area * reaction_area_factor
        
        # 计算扩散控制下的反应速率
        reaction_rate = k_m * (c_inf - c_s) * effective_area
        
        return reaction_rate
