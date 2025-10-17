"""
Comprehensive Safety Assessment Tool - 综合安全评估
Similar to ChemCrow's SafetySummary with enhanced features
"""

from typing import Dict, Any, List, Optional, Tuple
import re
from dataclasses import dataclass
from enum import Enum
import json
from pathlib import Path


class HazardClass(Enum):
    """GHS危害分类"""
    EXPLOSIVE = "Explosive"
    FLAMMABLE = "Flammable"
    OXIDIZING = "Oxidizing"
    CORROSIVE = "Corrosive"
    TOXIC = "Toxic"
    HARMFUL = "Harmful"
    ENVIRONMENTAL = "Environmental hazard"
    CARCINOGENIC = "Carcinogenic"
    MUTAGENIC = "Mutagenic"
    REPRODUCTIVE_TOXIC = "Reproductive toxicity"


@dataclass
class SafetyProfile:
    """安全概况数据结构"""
    ghs_classification: List[HazardClass]
    ghs_pictograms: List[str]
    signal_word: str  # "Danger" or "Warning"
    hazard_statements: List[str]
    precautionary_statements: List[str]
    ld50_oral: Optional[float]  # mg/kg
    ld50_dermal: Optional[float]
    lc50_inhalation: Optional[float]  # mg/L
    flash_point: Optional[float]  # °C
    autoignition_temp: Optional[float]  # °C
    explosive_limits: Optional[Tuple[float, float]]  # % in air
    environmental_hazards: Dict[str, Any]
    regulatory_status: Dict[str, Any]
    handling_recommendations: List[str]
    ppe_requirements: List[str]
    storage_conditions: List[str]
    disposal_methods: List[str]
    emergency_measures: Dict[str, str]


class SafetyAssessmentTool:
    """
    综合安全评估工具
    提供类似ChemCrow SafetySummary的功能，包括：
    - GHS分类
    - 爆炸性检查
    - 受控物质检查
    - 环境影响评估
    - 操作安全建议
    """
    
    def __init__(self):
        # 加载危险片段数据库
        self.explosive_fragments = self._load_explosive_fragments()
        self.toxic_fragments = self._load_toxic_fragments()
        self.controlled_substances = self._load_controlled_substances()
        self.environmental_hazards = self._load_environmental_hazards()
        
    def assess_safety(self, smiles: str, name: Optional[str] = None) -> Dict[str, Any]:
        """
        全面安全评估
        
        Args:
            smiles: SMILES字符串
            name: 化合物名称
            
        Returns:
            安全评估报告
        """
        report = {
            "compound": {
                "smiles": smiles,
                "name": name or "Unknown"
            },
            "safety_summary": "",
            "risk_level": "LOW",  # LOW, MODERATE, HIGH, EXTREME
            "ghs_classification": {},
            "physical_hazards": {},
            "health_hazards": {},
            "environmental_hazards": {},
            "regulatory_status": {},
            "handling_safety": {},
            "emergency_response": {},
            "recommendations": []
        }
        
        # 1. 物理危害评估
        report["physical_hazards"] = self._assess_physical_hazards(smiles)
        
        # 2. 健康危害评估
        report["health_hazards"] = self._assess_health_hazards(smiles)
        
        # 3. 环境危害评估
        report["environmental_hazards"] = self._assess_environmental_hazards(smiles)
        
        # 4. 法规状态检查
        report["regulatory_status"] = self._check_regulatory_status(smiles, name)
        
        # 5. GHS分类
        report["ghs_classification"] = self._classify_ghs(report)
        
        # 6. 操作安全建议
        report["handling_safety"] = self._generate_handling_recommendations(report)
        
        # 7. 应急响应
        report["emergency_response"] = self._generate_emergency_procedures(report)
        
        # 8. 计算总体风险等级
        report["risk_level"] = self._calculate_risk_level(report)
        
        # 9. 生成安全摘要
        report["safety_summary"] = self._generate_safety_summary(report)
        
        # 10. 生成建议
        report["recommendations"] = self._generate_recommendations(report)
        
        return report
        
    def _assess_physical_hazards(self, smiles: str) -> Dict[str, Any]:
        """评估物理危害"""
        hazards = {
            "explosive": False,
            "flammable": False,
            "oxidizing": False,
            "pyrophoric": False,
            "self_reactive": False,
            "details": []
        }
        
        # 检查爆炸性片段
        explosive_groups = self._check_explosive_fragments(smiles)
        if explosive_groups:
            hazards["explosive"] = True
            hazards["details"].append(f"Contains explosive groups: {', '.join(explosive_groups)}")
            
        # 检查易燃性（基于官能团）
        if self._check_flammability(smiles):
            hazards["flammable"] = True
            hazards["details"].append("Contains flammable functional groups")
            
        # 检查氧化性
        if self._check_oxidizing_potential(smiles):
            hazards["oxidizing"] = True
            hazards["details"].append("May act as an oxidizing agent")
            
        # 检查自反应性
        if self._check_self_reactivity(smiles):
            hazards["self_reactive"] = True
            hazards["details"].append("May undergo self-reactive decomposition")
            
        return hazards
        
    def _assess_health_hazards(self, smiles: str) -> Dict[str, Any]:
        """评估健康危害"""
        hazards = {
            "acute_toxicity": "Unknown",
            "skin_corrosion": False,
            "eye_damage": False,
            "respiratory_sensitization": False,
            "skin_sensitization": False,
            "mutagenicity": False,
            "carcinogenicity": False,
            "reproductive_toxicity": False,
            "specific_organ_toxicity": False,
            "aspiration_hazard": False,
            "details": []
        }
        
        # 检查毒性片段
        toxic_groups = self._check_toxic_fragments(smiles)
        if toxic_groups:
            hazards["acute_toxicity"] = "Likely toxic"
            hazards["details"].append(f"Contains toxic groups: {', '.join(toxic_groups)}")
            
        # 检查腐蚀性（强酸/强碱官能团）
        if self._check_corrosivity(smiles):
            hazards["skin_corrosion"] = True
            hazards["eye_damage"] = True
            hazards["details"].append("May cause severe skin burns and eye damage")
            
        # 检查致敏性
        if self._check_sensitization_potential(smiles):
            hazards["skin_sensitization"] = True
            hazards["details"].append("May cause allergic skin reaction")
            
        # 检查致癌/致突变性（基于已知结构警告）
        carcinogenic_alerts = self._check_carcinogenic_alerts(smiles)
        if carcinogenic_alerts:
            hazards["carcinogenicity"] = True
            hazards["mutagenicity"] = True
            hazards["details"].append(f"Structural alerts for carcinogenicity: {', '.join(carcinogenic_alerts)}")
            
        return hazards
        
    def _assess_environmental_hazards(self, smiles: str) -> Dict[str, Any]:
        """评估环境危害"""
        hazards = {
            "aquatic_toxicity": "Unknown",
            "persistence": "Unknown",
            "bioaccumulation": "Unknown",
            "ozone_depletion": False,
            "details": []
        }
        
        # 检查水生毒性（基于LogP和某些官能团）
        logp_estimate = self._estimate_logp(smiles)
        if logp_estimate > 4:
            hazards["bioaccumulation"] = "High potential"
            hazards["details"].append(f"High LogP ({logp_estimate:.1f}) suggests bioaccumulation potential")
            
        # 检查持久性（含卤素等）
        if self._check_persistence_indicators(smiles):
            hazards["persistence"] = "Likely persistent"
            hazards["details"].append("Contains persistent chemical groups")
            
        # 检查臭氧层破坏（含氟氯烃等）
        if self._check_ozone_depleting(smiles):
            hazards["ozone_depletion"] = True
            hazards["details"].append("May deplete ozone layer")
            
        return hazards
        
    def _check_regulatory_status(self, smiles: str, name: Optional[str]) -> Dict[str, Any]:
        """检查法规状态"""
        status = {
            "controlled_substance": False,
            "chemical_weapon": False,
            "precursor": False,
            "restricted": False,
            "regulations": [],
            "details": []
        }
        
        # 检查受控物质列表
        if name:
            name_lower = name.lower()
            for substance, info in self.controlled_substances.items():
                if substance.lower() in name_lower or name_lower in substance.lower():
                    status["controlled_substance"] = True
                    status["regulations"].append(info["regulation"])
                    status["details"].append(f"{substance}: {info['category']}")
                    
        # 检查化学武器相关
        if self._check_chemical_weapon_related(smiles):
            status["chemical_weapon"] = True
            status["regulations"].append("Chemical Weapons Convention")
            status["details"].append("Related to chemical weapons or precursors")
            
        # 检查前体化学品
        if self._check_precursor_chemicals(smiles):
            status["precursor"] = True
            status["regulations"].append("Precursor Control")
            status["details"].append("Listed as precursor chemical")
            
        return status
        
    def _classify_ghs(self, report: Dict[str, Any]) -> Dict[str, Any]:
        """生成GHS分类"""
        ghs = {
            "pictograms": [],
            "signal_word": "Warning",
            "hazard_statements": [],
            "precautionary_statements": [],
            "hazard_classes": []
        }
        
        # 基于物理危害
        if report["physical_hazards"]["explosive"]:
            ghs["pictograms"].append("GHS01 - Exploding bomb")
            ghs["signal_word"] = "Danger"
            ghs["hazard_statements"].append("H201 - Explosive; mass explosion hazard")
            ghs["hazard_classes"].append("Explosive")
            
        if report["physical_hazards"]["flammable"]:
            ghs["pictograms"].append("GHS02 - Flame")
            ghs["hazard_statements"].append("H225 - Highly flammable liquid and vapor")
            ghs["hazard_classes"].append("Flammable")
            
        # 基于健康危害
        if report["health_hazards"]["acute_toxicity"] == "Likely toxic":
            ghs["pictograms"].append("GHS06 - Skull and crossbones")
            ghs["signal_word"] = "Danger"
            ghs["hazard_statements"].append("H301 - Toxic if swallowed")
            ghs["hazard_classes"].append("Acute Toxicity")
            
        if report["health_hazards"]["skin_corrosion"]:
            ghs["pictograms"].append("GHS05 - Corrosion")
            ghs["signal_word"] = "Danger"
            ghs["hazard_statements"].append("H314 - Causes severe skin burns and eye damage")
            ghs["hazard_classes"].append("Corrosive")
            
        if report["health_hazards"]["carcinogenicity"]:
            ghs["pictograms"].append("GHS08 - Health hazard")
            ghs["signal_word"] = "Danger"
            ghs["hazard_statements"].append("H350 - May cause cancer")
            ghs["hazard_classes"].append("Carcinogenic")
            
        # 基于环境危害
        if report["environmental_hazards"]["aquatic_toxicity"] != "Unknown" or \
           report["environmental_hazards"]["bioaccumulation"] == "High potential":
            ghs["pictograms"].append("GHS09 - Environment")
            ghs["hazard_statements"].append("H410 - Very toxic to aquatic life with long lasting effects")
            ghs["hazard_classes"].append("Environmental Hazard")
            
        # 添加预防措施
        if ghs["hazard_classes"]:
            ghs["precautionary_statements"].extend([
                "P201 - Obtain special instructions before use",
                "P280 - Wear protective gloves/protective clothing/eye protection",
                "P301+P330+P331 - IF SWALLOWED: Rinse mouth. Do NOT induce vomiting",
                "P501 - Dispose of contents/container in accordance with regulations"
            ])
            
        return ghs
        
    def _generate_handling_recommendations(self, report: Dict[str, Any]) -> Dict[str, Any]:
        """生成操作安全建议"""
        handling = {
            "ppe": [],
            "engineering_controls": [],
            "work_practices": [],
            "storage": [],
            "incompatibilities": []
        }
        
        # PPE建议
        if report["risk_level"] in ["HIGH", "EXTREME"]:
            handling["ppe"].extend([
                "Full face shield",
                "Chemical-resistant gloves (nitrile or neoprene)",
                "Lab coat or chemical-resistant suit",
                "Closed-toe shoes",
                "Respiratory protection if aerosol generation possible"
            ])
        else:
            handling["ppe"].extend([
                "Safety glasses",
                "Lab coat",
                "Nitrile gloves",
                "Closed-toe shoes"
            ])
            
        # 工程控制
        if report["physical_hazards"]["explosive"] or report["health_hazards"]["acute_toxicity"] == "Likely toxic":
            handling["engineering_controls"].append("Use in fume hood")
            handling["engineering_controls"].append("Explosion-proof electrical equipment if explosive")
            
        # 工作实践
        handling["work_practices"].extend([
            "Avoid contact with skin and eyes",
            "Do not breathe dust/fume/gas/mist/vapors/spray",
            "Wash hands thoroughly after handling",
            "Do not eat, drink or smoke when using this product"
        ])
        
        # 存储条件
        if report["physical_hazards"]["flammable"]:
            handling["storage"].append("Store in cool, well-ventilated area away from heat sources")
            handling["storage"].append("Keep away from oxidizers")
            
        if report["physical_hazards"]["explosive"]:
            handling["storage"].append("Store in explosion-proof container")
            handling["storage"].append("Keep away from heat, sparks, and open flames")
            
        # 不相容性
        if report["physical_hazards"]["oxidizing"]:
            handling["incompatibilities"].append("Reducing agents")
            handling["incompatibilities"].append("Combustible materials")
            
        return handling
        
    def _generate_emergency_procedures(self, report: Dict[str, Any]) -> Dict[str, str]:
        """生成应急程序"""
        procedures = {
            "eye_contact": "Immediately flush eyes with plenty of water for at least 15 minutes. Get medical attention.",
            "skin_contact": "Remove contaminated clothing. Wash skin with soap and water. Get medical attention if irritation develops.",
            "inhalation": "Move to fresh air. If not breathing, give artificial respiration. Get medical attention.",
            "ingestion": "Do NOT induce vomiting. Rinse mouth with water. Never give anything by mouth to an unconscious person. Get medical attention immediately.",
            "fire": "Use appropriate media for surrounding fire. Cool containers with water spray.",
            "spill": "Evacuate area. Wear appropriate PPE. Absorb with inert material. Collect in appropriate container for disposal."
        }
        
        # 根据具体危害调整
        if report["health_hazards"]["skin_corrosion"]:
            procedures["skin_contact"] = "Immediately flush with large amounts of water for at least 30 minutes. Remove contaminated clothing while flushing. Get immediate medical attention."
            procedures["eye_contact"] = "Immediately flush with water for at least 30 minutes, lifting upper and lower eyelids. Get immediate medical attention."
            
        if report["physical_hazards"]["explosive"]:
            procedures["fire"] = "EVACUATE AREA IMMEDIATELY. Do not fight fire when fire reaches explosives. Cool containers from maximum distance."
            procedures["spill"] = "EVACUATE ALL PERSONNEL. Do not touch or walk through spilled material. Prevent entry into waterways, sewers, basements."
            
        return procedures
        
    def _calculate_risk_level(self, report: Dict[str, Any]) -> str:
        """计算总体风险等级"""
        score = 0
        
        # 物理危害评分
        if report["physical_hazards"]["explosive"]:
            score += 10
        if report["physical_hazards"]["flammable"]:
            score += 3
        if report["physical_hazards"]["oxidizing"]:
            score += 3
            
        # 健康危害评分
        if report["health_hazards"]["acute_toxicity"] == "Likely toxic":
            score += 5
        if report["health_hazards"]["carcinogenicity"]:
            score += 8
        if report["health_hazards"]["skin_corrosion"]:
            score += 4
            
        # 环境危害评分
        if report["environmental_hazards"]["bioaccumulation"] == "High potential":
            score += 3
        if report["environmental_hazards"]["persistence"] == "Likely persistent":
            score += 2
            
        # 法规评分
        if report["regulatory_status"]["controlled_substance"]:
            score += 5
        if report["regulatory_status"]["chemical_weapon"]:
            score += 10
            
        # 确定风险等级
        if score >= 15:
            return "EXTREME"
        elif score >= 10:
            return "HIGH"
        elif score >= 5:
            return "MODERATE"
        else:
            return "LOW"
            
    def _generate_safety_summary(self, report: Dict[str, Any]) -> str:
        """生成安全摘要"""
        summary_parts = []
        
        # 风险等级
        risk_emoji = {"LOW": "✅", "MODERATE": "⚠️", "HIGH": "🔴", "EXTREME": "☠️"}
        summary_parts.append(f"{risk_emoji[report['risk_level']]} Risk Level: {report['risk_level']}")
        
        # 主要危害
        hazards = []
        if report["physical_hazards"]["explosive"]:
            hazards.append("Explosive")
        if report["health_hazards"]["acute_toxicity"] == "Likely toxic":
            hazards.append("Toxic")
        if report["health_hazards"]["carcinogenicity"]:
            hazards.append("Carcinogenic")
        if report["regulatory_status"]["controlled_substance"]:
            hazards.append("Controlled Substance")
            
        if hazards:
            summary_parts.append(f"Main Hazards: {', '.join(hazards)}")
            
        # GHS信号词
        if report["ghs_classification"]["signal_word"]:
            summary_parts.append(f"GHS Signal Word: {report['ghs_classification']['signal_word']}")
            
        # 关键建议
        if report["risk_level"] in ["HIGH", "EXTREME"]:
            summary_parts.append("⚠️ SPECIAL PRECAUTIONS REQUIRED - Consult safety officer before use")
            
        return " | ".join(summary_parts)
        
    def _generate_recommendations(self, report: Dict[str, Any]) -> List[str]:
        """生成具体建议"""
        recommendations = []
        
        # 基于风险等级的通用建议
        if report["risk_level"] == "EXTREME":
            recommendations.append("🚫 Consider safer alternatives if possible")
            recommendations.append("📋 Require special training and approval before use")
            recommendations.append("🏥 Ensure medical surveillance for workers")
            
        elif report["risk_level"] == "HIGH":
            recommendations.append("⚠️ Implement strict controls and monitoring")
            recommendations.append("📚 Provide comprehensive safety training")
            recommendations.append("🧪 Use minimum quantities necessary")
            
        # 特定危害建议
        if report["physical_hazards"]["explosive"]:
            recommendations.append("💥 Store in explosion-proof facilities")
            recommendations.append("🌡️ Monitor temperature and avoid friction/impact")
            
        if report["health_hazards"]["carcinogenicity"]:
            recommendations.append("🔬 Designate as carcinogen area")
            recommendations.append("📝 Maintain exposure records")
            
        if report["environmental_hazards"]["bioaccumulation"] == "High potential":
            recommendations.append("♻️ Implement strict waste management")
            recommendations.append("🌊 Prevent release to environment")
            
        if report["regulatory_status"]["controlled_substance"]:
            recommendations.append("🔒 Secure storage with access control")
            recommendations.append("📊 Maintain detailed inventory records")
            
        return recommendations
        
    # ===== 辅助方法 =====
    
    def _load_explosive_fragments(self) -> List[str]:
        """加载爆炸性片段数据库"""
        return [
            "N(=O)O",  # Nitro
            "N=N=N",   # Azide
            "N#N+[O-]", # Diazo
            "C(=O)OO",  # Peroxide
            "N=N(=O)",  # Nitroso
            "C#N+[O-]", # Fulminates
            "[N+](=O)[O-]", # Nitrate ester
            "NN",       # Hydrazine derivatives
            "C(Cl)(Cl)(Cl)", # Polyhalogenated
            "[O-][O-]",  # Peroxide ion
        ]
        
    def _load_toxic_fragments(self) -> List[str]:
        """加载毒性片段数据库"""
        return [
            "C#N",      # Nitrile/Cyanide
            "[As]",     # Arsenic
            "[Hg]",     # Mercury
            "[Pb]",     # Lead
            "[Cd]",     # Cadmium
            "P(=S)",    # Thiophosphate
            "S(=O)(=O)F", # Sulfonyl fluoride
            "C(=O)Cl",  # Acyl chloride
            "N=C=O",    # Isocyanate
            "C=C(Cl)Cl", # Vinyl chloride
        ]
        
    def _load_controlled_substances(self) -> Dict[str, Dict[str, str]]:
        """加载受控物质数据库"""
        return {
            "Sarin": {"category": "Chemical Weapon", "regulation": "CWC Schedule 1"},
            "VX": {"category": "Chemical Weapon", "regulation": "CWC Schedule 1"},
            "Mustard gas": {"category": "Chemical Weapon", "regulation": "CWC Schedule 1"},
            "Ricin": {"category": "Toxin", "regulation": "CWC Schedule 1"},
            "Saxitoxin": {"category": "Toxin", "regulation": "CWC Schedule 1"},
            "Thiodiglycol": {"category": "Precursor", "regulation": "CWC Schedule 2"},
            "MDMA": {"category": "Controlled Drug", "regulation": "DEA Schedule I"},
            "Fentanyl": {"category": "Controlled Drug", "regulation": "DEA Schedule II"},
            "Ephedrine": {"category": "Precursor", "regulation": "DEA List I"},
            "Safrole": {"category": "Precursor", "regulation": "DEA List I"},
        }
        
    def _load_environmental_hazards(self) -> Dict[str, Any]:
        """加载环境危害数据"""
        return {
            "persistent_groups": ["C(F)", "C(Cl)", "Br", "[Si]"],
            "bioaccumulative_threshold_logp": 4.0,
            "ozone_depleting": ["C(F)(Cl)", "C(Br)"],
        }
        
    def _check_explosive_fragments(self, smiles: str) -> List[str]:
        """检查爆炸性片段"""
        found = []
        for fragment in self.explosive_fragments:
            if fragment in smiles:
                found.append(fragment)
        return found
        
    def _check_toxic_fragments(self, smiles: str) -> List[str]:
        """检查毒性片段"""
        found = []
        for fragment in self.toxic_fragments:
            if fragment in smiles:
                found.append(fragment)
        return found
        
    def _check_flammability(self, smiles: str) -> bool:
        """检查易燃性"""
        flammable_patterns = [
            r'CC',      # Hydrocarbons
            r'O[CH]',   # Ethers
            r'C=C',     # Alkenes
            r'C#C',     # Alkynes
            r'CO',      # Alcohols
        ]
        for pattern in flammable_patterns:
            if re.search(pattern, smiles):
                return True
        return False
        
    def _check_oxidizing_potential(self, smiles: str) -> bool:
        """检查氧化性"""
        oxidizing_patterns = [
            r'O[Cl]',   # Hypochlorites
            r'\[O-\]',  # Peroxides
            r'N\(=O\)=O', # Nitrates
            r'Cl\(=O\)', # Chlorates
        ]
        for pattern in oxidizing_patterns:
            if re.search(pattern, smiles):
                return True
        return False
        
    def _check_self_reactivity(self, smiles: str) -> bool:
        """检查自反应性"""
        # 检查是否同时含有氧化和还原基团
        has_oxidizing = self._check_oxidizing_potential(smiles)
        has_reducing = bool(re.search(r'[NH]|S', smiles))
        return has_oxidizing and has_reducing
        
    def _check_corrosivity(self, smiles: str) -> bool:
        """检查腐蚀性"""
        corrosive_patterns = [
            r'S\(=O\)\(=O\)O',  # Sulfonic acid
            r'P\(=O\)\(O\)O',   # Phosphoric acid
            r'C\(=O\)O',        # Carboxylic acid (strong ones)
            r'O[NH3+]',         # Strong bases
        ]
        for pattern in corrosive_patterns:
            if re.search(pattern, smiles):
                return True
        return False
        
    def _check_sensitization_potential(self, smiles: str) -> bool:
        """检查致敏性"""
        sensitizing_patterns = [
            r'N=C=O',    # Isocyanates
            r'C\(=O\)C=C', # α,β-unsaturated carbonyls
            r'C\(=O\)Cl', # Acyl chlorides
        ]
        for pattern in sensitizing_patterns:
            if re.search(pattern, smiles):
                return True
        return False
        
    def _check_carcinogenic_alerts(self, smiles: str) -> List[str]:
        """检查致癌性结构警告"""
        alerts = []
        carcinogenic_patterns = {
            "Aromatic amine": r'c[NH2]',
            "Nitrosamine": r'N\(N=O\)',
            "Polycyclic aromatic": r'c1ccc2c\(c1\)ccc1c2cccc1',
            "Alkylating agent": r'C\(Cl\)\(Cl\)',
            "Aromatic nitro": r'c\[N+\]\(=O\)\[O-\]',
        }
        for name, pattern in carcinogenic_patterns.items():
            if re.search(pattern, smiles):
                alerts.append(name)
        return alerts
        
    def _estimate_logp(self, smiles: str) -> float:
        """估算LogP（简化方法）"""
        # 这是一个非常简化的LogP估算
        carbon_count = smiles.count('C') + smiles.count('c')
        oxygen_count = smiles.count('O')
        nitrogen_count = smiles.count('N')
        
        # 简单的加权计算
        logp = carbon_count * 0.5 - oxygen_count * 0.5 - nitrogen_count * 0.3
        return min(max(logp, -3), 8)  # 限制在合理范围
        
    def _check_persistence_indicators(self, smiles: str) -> bool:
        """检查持久性指标"""
        for group in self.environmental_hazards["persistent_groups"]:
            if group in smiles:
                return True
        return False
        
    def _check_ozone_depleting(self, smiles: str) -> bool:
        """检查臭氧层破坏物质"""
        for group in self.environmental_hazards["ozone_depleting"]:
            if group in smiles:
                return True
        return False
        
    def _check_chemical_weapon_related(self, smiles: str) -> bool:
        """检查化学武器相关"""
        # 检查特定的化学武器相关结构
        cw_patterns = [
            r'P\(=O\)\(F\)',     # Nerve agents
            r'S\(CCCl\)\(CCCl\)', # Mustard gas
            r'C\(=O\)C\(Cl\)\(Cl\)Cl', # Phosgene related
        ]
        for pattern in cw_patterns:
            if re.search(pattern, smiles):
                return True
        return False
        
    def _check_precursor_chemicals(self, smiles: str) -> bool:
        """检查前体化学品"""
        precursor_patterns = [
            r'CC\(C\)NC',        # Ephedrine-like
            r'OCOc1ccc',         # Safrole-like
            r'P\(Cl\)\(Cl\)Cl',  # Phosphorus trichloride
        ]
        for pattern in precursor_patterns:
            if re.search(pattern, smiles):
                return True
        return False
