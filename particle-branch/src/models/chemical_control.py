import numpy as np

class ChemicalControl:
    def __init__(self):
        # 物理常数
        self.Ru = 8314.5  # 普适气体常数 (J/(kmol·K))
        
        # 反应相关参数
        self.h_reaction = -601.6e3  # Mg + CO2 -> MgO + C 的反应焓 (J/mol)
        self.M_mg = 24.305  # Mg的摩尔质量 (g/mol)
        self.MgO = 40.304
        self.M_co2 = 44.01  # CO2的摩尔质量 (g/mol)
        self.A_pre = 1.376e-4  # 指前因子 (s-1)
        self.E_a = 132.4e3  # 活化能 (J/mol)

        # 反应计量系数
        self.stoichiometry = {
            "reactants": {
                "Mg": 2,    # Mg的计量系数
                "CO2": 1    # CO2的计量系数
            },
            "products": {
                "MgO": 2,   # MgO的计量系数
                "C": 1      # C的计量系数
            }
        }

        # JANAF热力学数据
        # 格式：[T_low, T_high, T_common, [high_T_coeffs], [low_T_coeffs]]
        self.janaf_data = {
            "Mg": {
                "temp_range": [300.0, 3000.0, 1000.0],
                "high_T": [-2.403153, 7.37319e-2, -9.647178e-5, 6.288965e-8, -1.589644e-11, -4.851586e2, 7.0896e0],
                "low_T": [2.915891, 1.265919e-3, -4.297782e-7, 6.805673e-11, -3.985870e-15, -7.448851e2, -8.961338e-1]
            },
            "CO2": {
                "temp_range": [300.0, 5000.0, 1000.0],
                "high_T": [4.63659493, 2.74131991e-3, -9.95828531e-7, 1.60373011e-10, -9.16103468e-15, -4.89624199e4, -1.93534855],
                "low_T": [2.35677352, 8.98459677e-3, -7.12356269e-6, 2.45919022e-9, -1.43699548e-13, -4.83719697e4, 9.90105222]
            },
            "MgO": {
                "temp_range": [300.0, 5000.0, 1000.0],
                "high_T": [2.719002, 7.277419e-4, -1.455770e-7, 1.561599e-11, -6.373850e-16, -4.639325e4, 5.784988],
                "low_T": [2.507053, 1.726527e-3, -1.253360e-6, 4.473075e-10, -6.429272e-14, -4.632486e4, 6.571288]
            },
            "C": {
                "temp_range": [300.0, 5000.0, 1000.0],
                "high_T": [2.605583, -1.959343e-4, 1.067372e-7, -1.642394e-11, 8.187058e-16, 8.543678e4, 4.195308],
                "low_T": [2.554239, -3.215377e-4, 7.337087e-7, -7.322348e-10, 2.665214e-13, 8.547321e4, 4.531308]
            }
        }

    def get_coeffs(self, species, T):
        """获取对应温度下的JANAF系数"""
        if species not in self.janaf_data:
            raise ValueError(f"Species {species} not found in JANAF data")
        
        data = self.janaf_data[species]
        T_low, T_high, T_common = data["temp_range"]
        
        # 限制温度在有效范围内
        T = max(min(T, T_high), T_low)
        
        # 选择合适的温度范围系数
        if T < T_common:
            return data["low_T"]
        else:
            return data["high_T"]

    def calculate_Cp(self, species, T):
        """使用JANAF多项式计算比热容 (J/mol·K)"""
        a = self.get_coeffs(species, T)
        
        # 计算Cp/R
        Cp_R = a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4
        
        # 转换为J/mol·K
        return Cp_R * self.Ru

    def calculate_enthalpy(self, species, T):
        """计算特定温度下的摩尔焓 (J/mol)"""
        a = self.get_coeffs(species, T)
        
        # 计算H/RT
        H_RT = (a[0] + (a[1]/2)*T + (a[2]/3)*T**2 + (a[3]/4)*T**3 + (a[4]/5)*T**4 + a[5]/T)
        
        # 转换为J/mol
        return H_RT * self.Ru * T

    def calculate_entropy(self, species, T):
        """计算特定温度下的摩尔熵 (J/mol·K)"""
        a = self.get_coeffs(species, T)
        
        # 计算S/R
        S_R = (a[0] * np.log(T) + a[1]*T + (a[2]/2)*T**2 + (a[3]/3)*T**3 + (a[4]/4)*T**4 + a[6])
        
        # 转换为J/mol·K
        return S_R * self.Ru

    def calculate_gibbs(self, species, T):
        """计算特定温度下的摩尔吉布斯自由能 (J/mol)"""
        return self.calculate_enthalpy(species, T) - T * self.calculate_entropy(species, T)

    def calculate_reaction_enthalpy(self, T):
        """计算反应焓变，考虑化学计量数"""
        # 计算产物总焓
        H_products = sum(
            coeff * self.calculate_enthalpy(species, T)
            for species, coeff in self.stoichiometry["products"].items()
        )
        
        # 计算反应物总焓
        H_reactants = sum(
            coeff * self.calculate_enthalpy(species, T)
            for species, coeff in self.stoichiometry["reactants"].items()
        )
        
        # 返回反应焓变
        return H_products - H_reactants

    def calculate_reaction_gibbs(self, T):
        """计算反应吉布斯自由能变，考虑化学计量数"""
        # 计算产物总吉布斯自由能
        G_products = sum(
            coeff * self.calculate_gibbs(species, T)
            for species, coeff in self.stoichiometry["products"].items()
        )
        
        # 计算反应物总吉布斯自由能
        G_reactants = sum(
            coeff * self.calculate_gibbs(species, T)
            for species, coeff in self.stoichiometry["reactants"].items()
        )
        
        # 返回反应吉布斯自由能变
        return G_products - G_reactants

    def compute_chemical_controlled_rate(self, particle_temp, particle_diameter, 
                                      ambient_composition, ambient_pressure, reaction_area_factor):
        """计算化学动力学控制下的反应速率"""
        # 计算CO2分压
        p_co2 = ambient_pressure * ambient_composition.get("CO2", 0)
        
        # 计算反应速率常数（阿伦尼乌斯方程）
        k = self.A_pre * np.exp(-self.E_a / (self.Ru * particle_temp))
        
        # 计算反应面积
        surface_area = np.pi * particle_diameter**2
        effective_area = surface_area * reaction_area_factor
        
        # 计算反应速率（考虑CO2的化学计量数）
        reaction_rate = k * p_co2 * effective_area 
        
        return reaction_rate

    def calculate_reaction_rate(self, T, p_co2):
        """计算化学反应速率
        
        Args:
            T: 温度 (K)
            p_co2: CO2分压 (Pa)
            
        Returns:
            float: 反应速率 (kg/s)
        """
        # 更新有效反应面积系数
        self.update_reaction_area_factor()
        
        # 计算有效反应面积
        effective_area = self.calculate_effective_area()
        
        # 计算反应速率常数
        k = self.calculate_rate_constant(T)
        
        # 计算反应速率 (kg/s)
        rate = k * p_co2 * effective_area
        
        return rate
