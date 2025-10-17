"""
Chemical Price Lookup Tool - 化学品价格查询
Similar to ChemCrow's SMILES2Price
"""

from typing import Dict, Any, List, Optional, Tuple
import asyncio
import httpx
from datetime import datetime
from dataclasses import dataclass
from enum import Enum
import re
import json


class Supplier(Enum):
    """化学品供应商"""
    SIGMA_ALDRICH = "Sigma-Aldrich"
    TCI = "TCI Chemicals"
    ALFA_AESAR = "Alfa Aesar"
    ACROS = "Acros Organics"
    CHEMSPIDER = "ChemSpider"
    MCULE = "Mcule"
    MOLPORT = "MolPort"


@dataclass
class PriceInfo:
    """价格信息数据结构"""
    supplier: str
    catalog_number: str
    quantity: str
    unit: str
    price: float
    currency: str
    purity: str
    availability: str
    lead_time: Optional[str]
    url: Optional[str]
    last_updated: datetime
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "supplier": self.supplier,
            "catalog_number": self.catalog_number,
            "quantity": self.quantity,
            "unit": self.unit,
            "price": self.price,
            "currency": self.currency,
            "price_per_gram": self._calculate_price_per_gram(),
            "purity": self.purity,
            "availability": self.availability,
            "lead_time": self.lead_time,
            "url": self.url,
            "last_updated": self.last_updated.isoformat()
        }
        
    def _calculate_price_per_gram(self) -> float:
        """计算每克价格"""
        try:
            # 转换单位到克
            if self.unit.lower() == 'g':
                grams = float(re.findall(r'[\d.]+', self.quantity)[0])
            elif self.unit.lower() == 'mg':
                grams = float(re.findall(r'[\d.]+', self.quantity)[0]) / 1000
            elif self.unit.lower() == 'kg':
                grams = float(re.findall(r'[\d.]+', self.quantity)[0]) * 1000
            else:
                return 0
                
            return round(self.price / grams, 2)
        except:
            return 0


class PriceLookupTool:
    """
    化学品价格查询工具
    查询多个供应商的价格并提供购买建议
    """
    
    def __init__(self):
        self.timeout = httpx.Timeout(30.0)
        # 注：实际使用需要供应商API密钥
        self.api_keys = {
            Supplier.SIGMA_ALDRICH: None,
            Supplier.MCULE: None,
            Supplier.MOLPORT: None,
        }
        
    async def lookup_price(self,
                          smiles: Optional[str] = None,
                          name: Optional[str] = None,
                          cas: Optional[str] = None,
                          suppliers: Optional[List[Supplier]] = None) -> Dict[str, Any]:
        """
        查询化学品价格
        
        Args:
            smiles: SMILES字符串
            name: 化合物名称
            cas: CAS号
            suppliers: 指定供应商列表
            
        Returns:
            价格信息汇总
        """
        
        if not any([smiles, name, cas]):
            return {
                "error": "Please provide at least one identifier (SMILES, name, or CAS)"
            }
            
        # 默认查询所有主要供应商
        if suppliers is None:
            suppliers = [
                Supplier.SIGMA_ALDRICH,
                Supplier.TCI,
                Supplier.MCULE,
                Supplier.MOLPORT
            ]
            
        # 并行查询所有供应商
        tasks = []
        for supplier in suppliers:
            if supplier == Supplier.SIGMA_ALDRICH:
                tasks.append(self._query_sigma_aldrich(smiles, name, cas))
            elif supplier == Supplier.TCI:
                tasks.append(self._query_tci(smiles, name, cas))
            elif supplier == Supplier.MCULE:
                tasks.append(self._query_mcule(smiles, name, cas))
            elif supplier == Supplier.MOLPORT:
                tasks.append(self._query_molport(smiles, name, cas))
            else:
                tasks.append(self._query_generic_supplier(supplier, smiles, name, cas))
                
        results = await asyncio.gather(*tasks, return_exceptions=True)
        
        # 整理结果
        all_prices = []
        errors = []
        
        for i, result in enumerate(results):
            if isinstance(result, Exception):
                errors.append(f"{suppliers[i].value}: {str(result)}")
            elif result:
                all_prices.extend(result)
                
        # 分析价格数据
        analysis = self._analyze_prices(all_prices)
        
        # 生成购买建议
        recommendations = self._generate_recommendations(all_prices, analysis)
        
        return {
            "query": {
                "smiles": smiles,
                "name": name,
                "cas": cas
            },
            "prices": [p.to_dict() for p in all_prices],
            "total_suppliers": len(set(p.supplier for p in all_prices)),
            "total_offerings": len(all_prices),
            "analysis": analysis,
            "recommendations": recommendations,
            "errors": errors,
            "timestamp": datetime.now().isoformat()
        }
        
    async def _query_sigma_aldrich(self, smiles: Optional[str], 
                                  name: Optional[str], 
                                  cas: Optional[str]) -> List[PriceInfo]:
        """查询Sigma-Aldrich价格"""
        prices = []
        
        # 模拟API调用（实际需要真实API）
        # 这里返回示例数据
        if any([smiles, name, cas]):
            prices.append(PriceInfo(
                supplier=Supplier.SIGMA_ALDRICH.value,
                catalog_number="A12345",
                quantity="1",
                unit="g",
                price=125.00,
                currency="USD",
                purity="≥98%",
                availability="In Stock",
                lead_time=None,
                url="https://www.sigmaaldrich.com/catalog/product/aldrich/A12345",
                last_updated=datetime.now()
            ))
            
            prices.append(PriceInfo(
                supplier=Supplier.SIGMA_ALDRICH.value,
                catalog_number="A12345",
                quantity="5",
                unit="g",
                price=450.00,
                currency="USD",
                purity="≥98%",
                availability="In Stock",
                lead_time=None,
                url="https://www.sigmaaldrich.com/catalog/product/aldrich/A12345",
                last_updated=datetime.now()
            ))
            
        return prices
        
    async def _query_tci(self, smiles: Optional[str], 
                        name: Optional[str], 
                        cas: Optional[str]) -> List[PriceInfo]:
        """查询TCI价格"""
        prices = []
        
        # 模拟API调用
        if any([smiles, name, cas]):
            prices.append(PriceInfo(
                supplier=Supplier.TCI.value,
                catalog_number="T0123",
                quantity="1",
                unit="g",
                price=98.00,
                currency="USD",
                purity=">98.0%",
                availability="In Stock",
                lead_time="2-3 days",
                url="https://www.tcichemicals.com/US/en/p/T0123",
                last_updated=datetime.now()
            ))
            
        return prices
        
    async def _query_mcule(self, smiles: Optional[str], 
                          name: Optional[str], 
                          cas: Optional[str]) -> List[PriceInfo]:
        """查询Mcule价格（建筑块数据库）"""
        prices = []
        
        try:
            if smiles:
                # Mcule API示例（需要真实API密钥）
                async with httpx.AsyncClient(timeout=self.timeout) as client:
                    # 这是示例URL，实际需要正确的API端点
                    url = "https://mcule.com/api/v1/search"
                    headers = {"Authorization": f"Bearer {self.api_keys[Supplier.MCULE]}"}
                    params = {"smiles": smiles}
                    
                    # 模拟返回数据
                    prices.append(PriceInfo(
                        supplier=Supplier.MCULE.value,
                        catalog_number="MCULE-1234567890",
                        quantity="10",
                        unit="mg",
                        price=65.00,
                        currency="USD",
                        purity=">95%",
                        availability="Make-on-demand",
                        lead_time="4-6 weeks",
                        url="https://mcule.com/MCULE-1234567890",
                        last_updated=datetime.now()
                    ))
                    
        except Exception:
            pass
            
        return prices
        
    async def _query_molport(self, smiles: Optional[str], 
                            name: Optional[str], 
                            cas: Optional[str]) -> List[PriceInfo]:
        """查询MolPort价格（化合物市场）"""
        prices = []
        
        # 模拟API调用
        if any([smiles, name, cas]):
            prices.append(PriceInfo(
                supplier=Supplier.MOLPORT.value,
                catalog_number="MolPort-001-234-567",
                quantity="100",
                unit="mg",
                price=180.00,
                currency="EUR",
                purity=">90%",
                availability="2 in stock",
                lead_time="1-2 weeks",
                url="https://www.molport.com/shop/molecule-link/MolPort-001-234-567",
                last_updated=datetime.now()
            ))
            
        return prices
        
    async def _query_generic_supplier(self, supplier: Supplier,
                                     smiles: Optional[str], 
                                     name: Optional[str], 
                                     cas: Optional[str]) -> List[PriceInfo]:
        """通用供应商查询（占位）"""
        # 这里可以添加更多供应商的API集成
        return []
        
    def _analyze_prices(self, prices: List[PriceInfo]) -> Dict[str, Any]:
        """分析价格数据"""
        if not prices:
            return {"message": "No price data available"}
            
        analysis = {
            "price_range": {},
            "best_value": {},
            "availability_summary": {},
            "purity_levels": {},
            "quantity_options": []
        }
        
        # 转换为统一货币（USD）进行比较
        usd_prices = []
        for price in prices:
            usd_price = self._convert_to_usd(price.price, price.currency)
            price_per_gram = price._calculate_price_per_gram()
            if price_per_gram > 0:
                usd_prices.append({
                    "price": usd_price,
                    "price_per_gram": price_per_gram,
                    "supplier": price.supplier,
                    "quantity": price.quantity,
                    "unit": price.unit,
                    "purity": price.purity,
                    "availability": price.availability
                })
                
        if usd_prices:
            # 价格范围
            prices_per_gram = [p["price_per_gram"] for p in usd_prices]
            analysis["price_range"] = {
                "min": min(prices_per_gram),
                "max": max(prices_per_gram),
                "average": sum(prices_per_gram) / len(prices_per_gram),
                "currency": "USD/g"
            }
            
            # 最佳价值（考虑价格和纯度）
            best_value = min(usd_prices, key=lambda x: x["price_per_gram"])
            analysis["best_value"] = {
                "supplier": best_value["supplier"],
                "price_per_gram": best_value["price_per_gram"],
                "quantity": best_value["quantity"],
                "purity": best_value["purity"]
            }
            
        # 可用性汇总
        availability_count = {}
        for price in prices:
            status = "In Stock" if "stock" in price.availability.lower() else "On Demand"
            availability_count[status] = availability_count.get(status, 0) + 1
        analysis["availability_summary"] = availability_count
        
        # 纯度等级
        purity_levels = set()
        for price in prices:
            purity_levels.add(price.purity)
        analysis["purity_levels"] = sorted(list(purity_levels))
        
        # 数量选项
        quantity_options = []
        for price in prices:
            quantity_options.append(f"{price.quantity} {price.unit}")
        analysis["quantity_options"] = sorted(set(quantity_options))
        
        return analysis
        
    def _generate_recommendations(self, prices: List[PriceInfo], 
                                 analysis: Dict[str, Any]) -> List[str]:
        """生成购买建议"""
        recommendations = []
        
        if not prices:
            recommendations.append("❌ No commercial sources found. Consider custom synthesis.")
            return recommendations
            
        # 基于价格分析
        if "best_value" in analysis and analysis["best_value"]:
            best = analysis["best_value"]
            recommendations.append(
                f"💰 Best value: {best['supplier']} at ${best['price_per_gram']:.2f}/g "
                f"({best['quantity']}, {best['purity']})"
            )
            
        # 基于可用性
        if "availability_summary" in analysis:
            if analysis["availability_summary"].get("In Stock", 0) > 0:
                recommendations.append(
                    f"✅ {analysis['availability_summary']['In Stock']} suppliers have immediate availability"
                )
            else:
                recommendations.append(
                    "⏳ All suppliers require synthesis on demand (4-6 weeks typical)"
                )
                
        # 基于价格范围
        if "price_range" in analysis and analysis["price_range"]:
            price_range = analysis["price_range"]
            if price_range["max"] / price_range["min"] > 3:
                recommendations.append(
                    f"📊 Large price variation (${price_range['min']:.2f}-${price_range['max']:.2f}/g). "
                    f"Compare quality and lead times carefully."
                )
                
        # 批量购买建议
        if len(analysis.get("quantity_options", [])) > 3:
            recommendations.append(
                "📦 Multiple quantity options available. Consider bulk purchase for better unit price."
            )
            
        # 纯度建议
        purity_levels = analysis.get("purity_levels", [])
        if purity_levels:
            if any(">99" in p or "≥99" in p for p in purity_levels):
                recommendations.append(
                    "🔬 High purity options (>99%) available for analytical applications"
                )
            elif all("<95" in p or "≤95" in p for p in purity_levels):
                recommendations.append(
                    "⚠️ Only lower purity (<95%) available. May need purification for some applications."
                )
                
        # 供应商多样性
        if analysis.get("total_suppliers", 0) > 3:
            recommendations.append(
                "✅ Multiple suppliers available, reducing supply chain risk"
            )
        elif analysis.get("total_suppliers", 0) == 1:
            recommendations.append(
                "⚠️ Single supplier dependency. Consider alternative sources for critical projects."
            )
            
        return recommendations
        
    def _convert_to_usd(self, price: float, currency: str) -> float:
        """货币转换（简化，实际需要实时汇率）"""
        exchange_rates = {
            "USD": 1.0,
            "EUR": 1.1,
            "GBP": 1.25,
            "JPY": 0.007,
            "CNY": 0.14
        }
        return price * exchange_rates.get(currency, 1.0)


class ChemicalEconomicsAnalyzer:
    """
    化学经济性分析器
    评估合成vs购买的经济性
    """
    
    def __init__(self):
        self.price_tool = PriceLookupTool()
        
    async def analyze_make_vs_buy(self,
                                  target_smiles: str,
                                  required_quantity: float,
                                  quantity_unit: str = "g") -> Dict[str, Any]:
        """
        分析自制vs购买的经济性
        
        Args:
            target_smiles: 目标化合物SMILES
            required_quantity: 需要的数量
            quantity_unit: 数量单位
            
        Returns:
            经济性分析报告
        """
        
        # 查询购买价格
        buy_prices = await self.price_tool.lookup_price(smiles=target_smiles)
        
        # 估算合成成本（简化模型）
        synthesis_cost = self._estimate_synthesis_cost(
            target_smiles, 
            required_quantity,
            quantity_unit
        )
        
        # 转换数量到克
        quantity_in_grams = self._convert_to_grams(required_quantity, quantity_unit)
        
        # 计算购买成本
        buy_cost = self._calculate_buy_cost(buy_prices, quantity_in_grams)
        
        # 比较分析
        analysis = {
            "target": target_smiles,
            "required_quantity": f"{required_quantity} {quantity_unit}",
            "buy_option": buy_cost,
            "make_option": synthesis_cost,
            "recommendation": "",
            "factors_to_consider": []
        }
        
        # 生成建议
        if buy_cost["total_cost"] and synthesis_cost["total_cost"]:
            cost_ratio = buy_cost["total_cost"] / synthesis_cost["total_cost"]
            
            if cost_ratio > 2:
                analysis["recommendation"] = "🔬 MAKE: Synthesis is significantly cheaper"
            elif cost_ratio > 1.2:
                analysis["recommendation"] = "⚖️ MAKE: Synthesis offers modest cost savings"
            elif cost_ratio > 0.8:
                analysis["recommendation"] = "⚖️ NEUTRAL: Similar costs, consider other factors"
            else:
                analysis["recommendation"] = "💰 BUY: Purchasing is more economical"
                
        # 其他考虑因素
        analysis["factors_to_consider"] = [
            "⏰ Time: Synthesis requires 2-5 days vs immediate purchase",
            "🎯 Purity: Commercial sources often provide higher purity",
            "📊 Scale: Synthesis may be better for large quantities",
            "🔬 IP: Synthesis provides know-how and potential IP",
            "⚠️ Safety: Consider hazards of synthesis intermediates",
            "🌱 Sustainability: Evaluate environmental impact"
        ]
        
        return analysis
        
    def _estimate_synthesis_cost(self, smiles: str, quantity: float, unit: str) -> Dict[str, Any]:
        """估算合成成本（简化模型）"""
        # 这是一个非常简化的成本模型
        # 实际应该基于逆合成分析和原料价格
        
        quantity_in_grams = self._convert_to_grams(quantity, unit)
        
        # 基础成本因素
        raw_materials_cost = quantity_in_grams * 10  # $10/g原料成本（估算）
        labor_cost = 8 * 50  # 8小时 * $50/小时
        overhead = (raw_materials_cost + labor_cost) * 0.3  # 30%间接成本
        
        return {
            "total_cost": raw_materials_cost + labor_cost + overhead,
            "breakdown": {
                "raw_materials": raw_materials_cost,
                "labor": labor_cost,
                "overhead": overhead
            },
            "time_required": "2-5 days",
            "success_rate": "70-90%",
            "currency": "USD"
        }
        
    def _calculate_buy_cost(self, price_data: Dict[str, Any], quantity_grams: float) -> Dict[str, Any]:
        """计算购买成本"""
        if not price_data.get("prices"):
            return {"total_cost": None, "message": "No commercial source available"}
            
        # 找最佳价格
        best_price_per_gram = price_data["analysis"].get("price_range", {}).get("min")
        
        if best_price_per_gram:
            total_cost = best_price_per_gram * quantity_grams
            
            return {
                "total_cost": total_cost,
                "unit_price": best_price_per_gram,
                "supplier": price_data["analysis"]["best_value"]["supplier"],
                "lead_time": "1-4 weeks typical",
                "currency": "USD"
            }
            
        return {"total_cost": None, "message": "Unable to calculate cost"}
        
    def _convert_to_grams(self, quantity: float, unit: str) -> float:
        """转换到克"""
        conversions = {
            "g": 1,
            "mg": 0.001,
            "kg": 1000,
            "lb": 453.592,
            "oz": 28.3495
        }
        return quantity * conversions.get(unit.lower(), 1)
