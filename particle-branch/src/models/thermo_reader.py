"""
读取热力学数据文件并计算热力学性质
"""

import os
import numpy as np

class ThermoReader:
    def __init__(self, thermo_file):
        """初始化热力学数据读取器
        
        Args:
            thermo_file: thermo.dat文件的路径
        """
        self.thermo_file = thermo_file
        self.species_data = {}
        self.read_thermo_file()
    
    def read_thermo_file(self):
        """读取thermo.dat文件并解析数据"""
        with open(self.thermo_file, 'r') as f:
            lines = f.readlines()
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line or line.startswith('!'): 
                i += 1
                continue
                
            # 读取物种信息行
            species_name = line[:18].strip()
            phase = line[44] if len(line) > 44 else 'G'  # 默认为气体
            
            try:
                # 读取三行系数
                i += 1
                cp_coeffs = [float(x) for x in lines[i].strip().split()]
                i += 1
                h_coeffs = [float(x) for x in lines[i].strip().split()]
                i += 1
                s_coeffs = [float(x) for x in lines[i].strip().split()]
                
                # 存储物种数据
                self.species_data[species_name] = {
                    'phase': phase,
                    'cp_coeffs': np.array(cp_coeffs),
                    'h_coeffs': np.array(h_coeffs),
                    's_coeffs': np.array(s_coeffs)
                }
            except (ValueError, IndexError) as e:
                print(f"警告：解析{species_name}的热力学数据时出错：{str(e)}")
            
            i += 1
    
    def calculate_Cp(self, species, T):
        """计算定压比热容
        
        Args:
            species: 物种名称
            T: 温度 (K)
            
        Returns:
            Cp: 定压比热容 (J/kg·K)
            
        Raises:
            ValueError: 如果物种不存在于热力学数据库中
        """
        if species not in self.species_data:
            raise ValueError(f"缺少物种 {species} 的热力学系数，无法计算比热容")
            
        a = self.species_data[species]['cp_coeffs']
        T_powers = np.array([1, T, T**2, T**3, T**4])
        return np.sum(a * T_powers)*self.Ru
    
    def calculate_H(self, species, T):
        """计算焓值
        
        Args:
            species: 物种名称
            T: 温度 (K)
            
        Returns:
            H: 焓值 (J/kg)
            
        Raises:
            ValueError: 如果物种不存在于热力学数据库中
        """
        if species not in self.species_data:
            raise ValueError(f"缺少物种 {species} 的热力学系数，无法计算焓值")
            
        a = self.species_data[species]['h_coeffs']
        T_powers = np.array([1, 1/2*T, 1/3*T**2, 1/4*T**3, 1/5*T**4])
        return np.sum(a * T_powers)*self.Ru*T
    
    def calculate_S(self, species, T):
        """计算熵
        
        Args:
            species: 物种名称
            T: 温度 (K)
            
        Returns:
            S: 熵 (J/kg·K)
            
        Raises:
            ValueError: 如果物种不存在于热力学数据库中
        """
        if species not in self.species_data:
            raise ValueError(f"缺少物种 {species} 的热力学系数，无法计算熵")
            
        a = self.species_data[species]['s_coeffs']
        T_powers = np.array([1, T, T**2, T**3, T**4])
        return np.sum(a * T_powers)
    
    def has_species(self, species):
        """检查是否存在指定物种的热力学数据
        
        Args:
            species: 物种名称
            
        Returns:
            bool: 是否存在该物种的数据
        """
        return species in self.species_data 