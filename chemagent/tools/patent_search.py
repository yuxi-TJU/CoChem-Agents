"""
Patent Search Tool - 专利检查功能
Searches USPTO, Google Patents, and other patent databases
"""

import asyncio
from typing import Dict, Any, List, Optional
import httpx
from datetime import datetime
import json
import re


class PatentSearchTool:
    """
    专利搜索工具 - 类似ChemCrow的PatentCheck
    检查分子是否已被专利保护
    """
    
    def __init__(self):
        self.uspto_api = "https://developer.uspto.gov/ibd-api/v1/patent/application"
        self.google_patents_api = "https://patents.google.com/xhr/query"
        self.timeout = httpx.Timeout(30.0)
        
    async def check_patent(self, 
                          smiles: Optional[str] = None,
                          name: Optional[str] = None,
                          cas: Optional[str] = None) -> Dict[str, Any]:
        """
        检查化合物是否有专利保护
        
        Args:
            smiles: SMILES字符串
            name: 化合物名称
            cas: CAS号
            
        Returns:
            专利信息字典
        """
        results = {
            "has_patents": False,
            "patent_count": 0,
            "patents": [],
            "search_terms": [],
            "warnings": []
        }
        
        # 构建搜索词
        search_terms = []
        if smiles:
            # 简化SMILES用于搜索（去除立体化学等）
            simplified = self._simplify_smiles_for_search(smiles)
            search_terms.append(simplified)
            
        if name:
            search_terms.append(name)
            # 添加常见变体
            search_terms.extend(self._generate_name_variants(name))
            
        if cas:
            search_terms.append(cas)
            # CAS格式：12345-67-8
            search_terms.append(cas.replace("-", ""))
            
        results["search_terms"] = search_terms
        
        # 并行搜索多个数据库
        tasks = []
        for term in search_terms[:3]:  # 限制搜索数量
            tasks.append(self._search_uspto(term))
            tasks.append(self._search_google_patents(term))
            
        search_results = await asyncio.gather(*tasks, return_exceptions=True)
        
        # 合并结果
        all_patents = []
        for result in search_results:
            if isinstance(result, Exception):
                results["warnings"].append(str(result))
            elif result:
                all_patents.extend(result)
                
        # 去重
        seen = set()
        unique_patents = []
        for patent in all_patents:
            patent_id = patent.get("patent_number", patent.get("title", ""))
            if patent_id not in seen:
                seen.add(patent_id)
                unique_patents.append(patent)
                
        results["patents"] = unique_patents[:20]  # 限制返回数量
        results["patent_count"] = len(unique_patents)
        results["has_patents"] = len(unique_patents) > 0
        
        # 分析专利状态
        if results["has_patents"]:
            results["analysis"] = self._analyze_patents(unique_patents)
            
        return results
        
    async def _search_uspto(self, query: str) -> List[Dict[str, Any]]:
        """搜索USPTO数据库"""
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                # USPTO API需要特定格式
                params = {
                    "searchText": query,
                    "start": 0,
                    "rows": 10
                }
                
                response = await client.get(
                    self.uspto_api,
                    params=params
                )
                
                if response.status_code == 200:
                    data = response.json()
                    patents = []
                    
                    for doc in data.get("response", {}).get("docs", []):
                        patents.append({
                            "source": "USPTO",
                            "patent_number": doc.get("patentNumber"),
                            "title": doc.get("inventionTitle"),
                            "assignee": doc.get("assigneeEntityName"),
                            "filing_date": doc.get("filingDate"),
                            "grant_date": doc.get("grantDate"),
                            "abstract": doc.get("abstractText", "")[:500],
                            "status": doc.get("applicationStatus"),
                            "url": f"https://patents.uspto.gov/patent/{doc.get('patentNumber')}"
                        })
                        
                    return patents
                    
        except Exception as e:
            # 返回空列表而不是抛出异常
            return []
            
        return []
        
    async def _search_google_patents(self, query: str) -> List[Dict[str, Any]]:
        """搜索Google Patents"""
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                # Google Patents使用不同的API格式
                params = {
                    "q": query,
                    "limit": 10,
                    "lang": "en"
                }
                
                response = await client.get(
                    self.google_patents_api,
                    params=params
                )
                
                if response.status_code == 200:
                    data = response.json()
                    patents = []
                    
                    for result in data.get("results", []):
                        patents.append({
                            "source": "Google Patents",
                            "patent_number": result.get("patent_number"),
                            "title": result.get("title"),
                            "assignee": result.get("assignee"),
                            "filing_date": result.get("filing_date"),
                            "publication_date": result.get("publication_date"),
                            "abstract": result.get("snippet", "")[:500],
                            "url": result.get("url")
                        })
                        
                    return patents
                    
        except Exception:
            return []
            
        return []
        
    def _simplify_smiles_for_search(self, smiles: str) -> str:
        """简化SMILES用于专利搜索"""
        # 移除立体化学信息
        simplified = re.sub(r'[@\\\/]', '', smiles)
        # 移除氢原子
        simplified = re.sub(r'\[H\]', '', simplified)
        # 简化芳香环
        simplified = simplified.replace('c', 'C')
        simplified = simplified.replace('n', 'N')
        simplified = simplified.replace('o', 'O')
        simplified = simplified.replace('s', 'S')
        
        return simplified
        
    def _generate_name_variants(self, name: str) -> List[str]:
        """生成化合物名称的变体用于搜索"""
        variants = []
        
        # 小写版本
        variants.append(name.lower())
        
        # 去除常见后缀
        for suffix in ['acid', 'ester', 'salt', 'hydrochloride', 'sulfate']:
            if name.lower().endswith(suffix):
                base_name = name[:-len(suffix)].strip()
                variants.append(base_name)
                
        # 处理连字符
        if '-' in name:
            variants.append(name.replace('-', ' '))
            variants.append(name.replace('-', ''))
            
        return variants[:3]  # 限制变体数量
        
    def _analyze_patents(self, patents: List[Dict[str, Any]]) -> Dict[str, Any]:
        """分析专利信息"""
        analysis = {
            "active_patents": 0,
            "expired_patents": 0,
            "pending_applications": 0,
            "main_assignees": [],
            "earliest_filing": None,
            "latest_filing": None,
            "geographic_coverage": set(),
            "recommendation": ""
        }
        
        assignee_count = {}
        
        for patent in patents:
            # 统计状态
            status = patent.get("status", "").lower()
            if "granted" in status or "active" in status:
                analysis["active_patents"] += 1
            elif "expired" in status:
                analysis["expired_patents"] += 1
            elif "pending" in status or "application" in status:
                analysis["pending_applications"] += 1
                
            # 统计申请人
            assignee = patent.get("assignee")
            if assignee:
                assignee_count[assignee] = assignee_count.get(assignee, 0) + 1
                
            # 记录日期
            filing_date = patent.get("filing_date")
            if filing_date:
                if not analysis["earliest_filing"] or filing_date < analysis["earliest_filing"]:
                    analysis["earliest_filing"] = filing_date
                if not analysis["latest_filing"] or filing_date > analysis["latest_filing"]:
                    analysis["latest_filing"] = filing_date
                    
        # 主要申请人
        analysis["main_assignees"] = sorted(
            assignee_count.items(), 
            key=lambda x: x[1], 
            reverse=True
        )[:3]
        
        # 生成建议
        if analysis["active_patents"] > 0:
            analysis["recommendation"] = "⚠️ 该化合物可能受专利保护，建议进行详细的自由实施(FTO)分析"
        elif analysis["expired_patents"] > 0:
            analysis["recommendation"] = "✅ 发现已过期专利，该化合物可能可以自由使用"
        elif analysis["pending_applications"] > 0:
            analysis["recommendation"] = "⏳ 存在待审专利申请，需要持续监控"
        else:
            analysis["recommendation"] = "✅ 未发现相关专利，但建议进行更全面的搜索"
            
        return analysis


class CASLookupTool:
    """
    CAS号查询工具 - 类似ChemCrow的Name2CAS
    """
    
    def __init__(self):
        self.cas_api = "https://commonchemistry.cas.org/api/search"
        
    async def name_to_cas(self, name: str) -> Dict[str, Any]:
        """
        将化合物名称转换为CAS号
        
        Args:
            name: 化合物名称（通用名、IUPAC名等）
            
        Returns:
            CAS信息字典
        """
        try:
            async with httpx.AsyncClient() as client:
                response = await client.get(
                    self.cas_api,
                    params={"q": name}
                )
                
                if response.status_code == 200:
                    data = response.json()
                    
                    if data.get("results"):
                        result = data["results"][0]
                        return {
                            "success": True,
                            "cas_number": result.get("rn"),
                            "name": result.get("name"),
                            "synonyms": result.get("synonyms", []),
                            "molecular_formula": result.get("molecularFormula"),
                            "molecular_weight": result.get("molecularMass"),
                            "smiles": result.get("smile"),
                            "inchi": result.get("inchi"),
                            "inchi_key": result.get("inchiKey")
                        }
                        
        except Exception as e:
            return {
                "success": False,
                "error": str(e),
                "cas_number": None
            }
            
        return {
            "success": False,
            "error": "CAS number not found",
            "cas_number": None
        }
        
    async def cas_to_name(self, cas: str) -> Dict[str, Any]:
        """
        将CAS号转换为化合物名称
        
        Args:
            cas: CAS注册号
            
        Returns:
            化合物信息字典
        """
        # 验证CAS格式
        if not re.match(r'^\d{2,7}-\d{2}-\d$', cas):
            return {
                "success": False,
                "error": "Invalid CAS number format",
                "name": None
            }
            
        try:
            async with httpx.AsyncClient() as client:
                # 直接查询CAS号
                response = await client.get(
                    f"https://commonchemistry.cas.org/api/detail",
                    params={"cas_rn": cas}
                )
                
                if response.status_code == 200:
                    data = response.json()
                    return {
                        "success": True,
                        "cas_number": cas,
                        "name": data.get("name"),
                        "synonyms": data.get("synonyms", []),
                        "molecular_formula": data.get("molecularFormula"),
                        "molecular_weight": data.get("molecularMass"),
                        "smiles": data.get("smile"),
                        "inchi": data.get("inchi")
                    }
                    
        except Exception as e:
            return {
                "success": False,
                "error": str(e),
                "name": None
            }
            
        return {
            "success": False,
            "error": "CAS information not found",
            "name": None
        }
