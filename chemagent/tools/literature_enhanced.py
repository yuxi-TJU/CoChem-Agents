"""
Enhanced Literature Search Tool - 增强版文献搜索
Similar to ChemCrow's LitSearch with vector search capabilities
"""

from typing import Dict, Any, List, Optional, Tuple
import asyncio
import httpx
from datetime import datetime
import xml.etree.ElementTree as ET
import json
import re
from urllib.parse import quote
import hashlib
from dataclasses import dataclass
from enum import Enum


class SearchSource(Enum):
    """文献搜索源"""
    PUBMED = "pubmed"
    CHEMRXIV = "chemrxiv"
    ARXIV = "arxiv"
    CROSSREF = "crossref"
    SEMANTIC_SCHOLAR = "semantic_scholar"
    GOOGLE_SCHOLAR = "google_scholar"


@dataclass
class Paper:
    """文献数据结构"""
    title: str
    authors: List[str]
    abstract: str
    doi: Optional[str]
    pmid: Optional[str]
    arxiv_id: Optional[str]
    year: int
    journal: Optional[str]
    url: str
    source: SearchSource
    citations: Optional[int] = None
    relevance_score: float = 0.0
    keywords: List[str] = None
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "title": self.title,
            "authors": self.authors,
            "abstract": self.abstract[:500] + "..." if len(self.abstract) > 500 else self.abstract,
            "doi": self.doi,
            "pmid": self.pmid,
            "year": self.year,
            "journal": self.journal,
            "url": self.url,
            "source": self.source.value,
            "citations": self.citations,
            "relevance_score": self.relevance_score,
            "keywords": self.keywords
        }


class EnhancedLiteratureSearch:
    """
    增强版文献搜索工具
    支持多数据源、相关性排序、引用分析
    """
    
    def __init__(self):
        self.timeout = httpx.Timeout(30.0)
        self.cache = {}  # 简单缓存
        
        # API endpoints
        self.endpoints = {
            SearchSource.PUBMED: "https://eutils.ncbi.nlm.nih.gov/entrez/eutils",
            SearchSource.CHEMRXIV: "https://chemrxiv.org/engage/chemrxiv/public-api/v1",
            SearchSource.ARXIV: "http://export.arxiv.org/api/query",
            SearchSource.CROSSREF: "https://api.crossref.org/works",
            SearchSource.SEMANTIC_SCHOLAR: "https://api.semanticscholar.org/v1/paper",
        }
        
    async def search(self,
                    query: str,
                    sources: List[SearchSource] = None,
                    max_results: int = 20,
                    year_from: Optional[int] = None,
                    year_to: Optional[int] = None,
                    sort_by: str = "relevance") -> Dict[str, Any]:
        """
        综合文献搜索
        
        Args:
            query: 搜索查询
            sources: 搜索源列表
            max_results: 最大结果数
            year_from: 起始年份
            year_to: 结束年份
            sort_by: 排序方式 (relevance, year, citations)
            
        Returns:
            搜索结果字典
        """
        
        # 默认搜索所有主要源
        if sources is None:
            sources = [
                SearchSource.PUBMED,
                SearchSource.CHEMRXIV,
                SearchSource.ARXIV,
                SearchSource.CROSSREF
            ]
            
        # 检查缓存
        cache_key = self._get_cache_key(query, sources, year_from, year_to)
        if cache_key in self.cache:
            cached_result = self.cache[cache_key]
            if (datetime.now() - cached_result["timestamp"]).seconds < 3600:
                return cached_result["data"]
                
        # 并行搜索所有源
        tasks = []
        for source in sources:
            if source == SearchSource.PUBMED:
                tasks.append(self._search_pubmed(query, max_results, year_from, year_to))
            elif source == SearchSource.CHEMRXIV:
                tasks.append(self._search_chemrxiv(query, max_results))
            elif source == SearchSource.ARXIV:
                tasks.append(self._search_arxiv(query, max_results))
            elif source == SearchSource.CROSSREF:
                tasks.append(self._search_crossref(query, max_results, year_from, year_to))
            elif source == SearchSource.SEMANTIC_SCHOLAR:
                tasks.append(self._search_semantic_scholar(query, max_results))
                
        results = await asyncio.gather(*tasks, return_exceptions=True)
        
        # 合并结果
        all_papers = []
        errors = []
        
        for i, result in enumerate(results):
            if isinstance(result, Exception):
                errors.append(f"{sources[i].value}: {str(result)}")
            else:
                all_papers.extend(result)
                
        # 去重（基于标题相似度）
        unique_papers = self._deduplicate_papers(all_papers)
        
        # 计算相关性分数
        for paper in unique_papers:
            paper.relevance_score = self._calculate_relevance(paper, query)
            
        # 排序
        if sort_by == "relevance":
            unique_papers.sort(key=lambda p: p.relevance_score, reverse=True)
        elif sort_by == "year":
            unique_papers.sort(key=lambda p: p.year, reverse=True)
        elif sort_by == "citations":
            unique_papers.sort(key=lambda p: p.citations or 0, reverse=True)
            
        # 限制结果数量
        unique_papers = unique_papers[:max_results]
        
        # 构建结果
        result = {
            "query": query,
            "total_found": len(unique_papers),
            "sources_searched": [s.value for s in sources],
            "papers": [p.to_dict() for p in unique_papers],
            "errors": errors,
            "search_time": datetime.now().isoformat(),
            "statistics": self._calculate_statistics(unique_papers)
        }
        
        # 缓存结果
        self.cache[cache_key] = {
            "timestamp": datetime.now(),
            "data": result
        }
        
        return result
        
    async def _search_pubmed(self, query: str, max_results: int, 
                           year_from: Optional[int], year_to: Optional[int]) -> List[Paper]:
        """搜索PubMed数据库"""
        papers = []
        
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                # 构建查询
                search_query = query
                if year_from or year_to:
                    year_from = year_from or 1900
                    year_to = year_to or datetime.now().year
                    search_query += f" AND {year_from}:{year_to}[dp]"
                    
                # 搜索获取ID列表
                search_url = f"{self.endpoints[SearchSource.PUBMED]}/esearch.fcgi"
                search_params = {
                    "db": "pubmed",
                    "term": search_query,
                    "retmax": max_results,
                    "retmode": "json"
                }
                
                search_response = await client.get(search_url, params=search_params)
                search_data = search_response.json()
                
                id_list = search_data.get("esearchresult", {}).get("idlist", [])
                
                if id_list:
                    # 获取详细信息
                    fetch_url = f"{self.endpoints[SearchSource.PUBMED]}/efetch.fcgi"
                    fetch_params = {
                        "db": "pubmed",
                        "id": ",".join(id_list),
                        "retmode": "xml"
                    }
                    
                    fetch_response = await client.get(fetch_url, params=fetch_params)
                    
                    # 解析XML
                    root = ET.fromstring(fetch_response.text)
                    
                    for article in root.findall(".//PubmedArticle"):
                        try:
                            # 提取信息
                            pmid = article.find(".//PMID").text
                            title = article.find(".//ArticleTitle").text
                            
                            # 作者
                            authors = []
                            for author in article.findall(".//Author"):
                                last_name = author.find("LastName")
                                fore_name = author.find("ForeName")
                                if last_name is not None and fore_name is not None:
                                    authors.append(f"{fore_name.text} {last_name.text}")
                                    
                            # 摘要
                            abstract_elem = article.find(".//AbstractText")
                            abstract = abstract_elem.text if abstract_elem is not None else ""
                            
                            # 期刊和年份
                            journal_elem = article.find(".//Journal/Title")
                            journal = journal_elem.text if journal_elem is not None else ""
                            
                            year_elem = article.find(".//PubDate/Year")
                            year = int(year_elem.text) if year_elem is not None else datetime.now().year
                            
                            # DOI
                            doi_elem = article.find(".//ArticleId[@IdType='doi']")
                            doi = doi_elem.text if doi_elem is not None else None
                            
                            # 关键词
                            keywords = [kw.text for kw in article.findall(".//Keyword")]
                            
                            papers.append(Paper(
                                title=title,
                                authors=authors[:5],  # 限制作者数量
                                abstract=abstract,
                                doi=doi,
                                pmid=pmid,
                                arxiv_id=None,
                                year=year,
                                journal=journal,
                                url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                                source=SearchSource.PUBMED,
                                keywords=keywords
                            ))
                            
                        except Exception:
                            continue
                            
        except Exception as e:
            # 记录错误但不中断
            pass
            
        return papers
        
    async def _search_chemrxiv(self, query: str, max_results: int) -> List[Paper]:
        """搜索ChemRxiv预印本"""
        papers = []
        
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                # ChemRxiv API
                url = f"{self.endpoints[SearchSource.CHEMRXIV]}/items"
                params = {
                    "term": query,
                    "limit": max_results,
                    "skip": 0
                }
                
                response = await client.get(url, params=params)
                
                if response.status_code == 200:
                    data = response.json()
                    
                    for item in data.get("itemHits", []):
                        try:
                            paper_data = item.get("item", {})
                            
                            papers.append(Paper(
                                title=paper_data.get("title", ""),
                                authors=[a.get("name", "") for a in paper_data.get("authors", [])[:5]],
                                abstract=paper_data.get("abstract", ""),
                                doi=paper_data.get("doi"),
                                pmid=None,
                                arxiv_id=None,
                                year=int(paper_data.get("publishedDate", "2024")[:4]),
                                journal="ChemRxiv (Preprint)",
                                url=paper_data.get("asset", {}).get("original", {}).get("url", ""),
                                source=SearchSource.CHEMRXIV,
                                keywords=paper_data.get("keywords", [])
                            ))
                            
                        except Exception:
                            continue
                            
        except Exception:
            pass
            
        return papers
        
    async def _search_arxiv(self, query: str, max_results: int) -> List[Paper]:
        """搜索arXiv预印本"""
        papers = []
        
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                # arXiv API
                params = {
                    "search_query": f"all:{query} OR cat:physics.chem-ph",
                    "start": 0,
                    "max_results": max_results
                }
                
                response = await client.get(self.endpoints[SearchSource.ARXIV], params=params)
                
                if response.status_code == 200:
                    # 解析Atom feed
                    root = ET.fromstring(response.text)
                    
                    # 定义命名空间
                    ns = {
                        'atom': 'http://www.w3.org/2005/Atom',
                        'arxiv': 'http://arxiv.org/schemas/atom'
                    }
                    
                    for entry in root.findall('atom:entry', ns):
                        try:
                            title = entry.find('atom:title', ns).text.strip()
                            
                            authors = []
                            for author in entry.findall('atom:author', ns):
                                name = author.find('atom:name', ns).text
                                authors.append(name)
                                
                            abstract = entry.find('atom:summary', ns).text.strip()
                            
                            # 提取arXiv ID
                            arxiv_id = entry.find('atom:id', ns).text
                            arxiv_id = arxiv_id.split('/')[-1]
                            
                            # 发布日期
                            published = entry.find('atom:published', ns).text
                            year = int(published[:4])
                            
                            # DOI（如果有）
                            doi_elem = entry.find('arxiv:doi', ns)
                            doi = doi_elem.text if doi_elem is not None else None
                            
                            papers.append(Paper(
                                title=title,
                                authors=authors[:5],
                                abstract=abstract,
                                doi=doi,
                                pmid=None,
                                arxiv_id=arxiv_id,
                                year=year,
                                journal="arXiv (Preprint)",
                                url=f"https://arxiv.org/abs/{arxiv_id}",
                                source=SearchSource.ARXIV
                            ))
                            
                        except Exception:
                            continue
                            
        except Exception:
            pass
            
        return papers
        
    async def _search_crossref(self, query: str, max_results: int,
                              year_from: Optional[int], year_to: Optional[int]) -> List[Paper]:
        """搜索CrossRef数据库"""
        papers = []
        
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                # CrossRef API
                params = {
                    "query": query,
                    "rows": max_results
                }
                
                # 添加年份过滤
                filters = []
                if year_from:
                    filters.append(f"from-pub-date:{year_from}")
                if year_to:
                    filters.append(f"until-pub-date:{year_to}")
                    
                if filters:
                    params["filter"] = ",".join(filters)
                    
                response = await client.get(
                    self.endpoints[SearchSource.CROSSREF],
                    params=params
                )
                
                if response.status_code == 200:
                    data = response.json()
                    
                    for item in data.get("message", {}).get("items", []):
                        try:
                            # 提取作者
                            authors = []
                            for author in item.get("author", [])[:5]:
                                given = author.get("given", "")
                                family = author.get("family", "")
                                authors.append(f"{given} {family}".strip())
                                
                            # 提取年份
                            date_parts = item.get("published-print", {}).get("date-parts", [[]])
                            if not date_parts:
                                date_parts = item.get("published-online", {}).get("date-parts", [[]])
                            year = date_parts[0][0] if date_parts and date_parts[0] else datetime.now().year
                            
                            papers.append(Paper(
                                title=item.get("title", [""])[0],
                                authors=authors,
                                abstract=item.get("abstract", ""),
                                doi=item.get("DOI"),
                                pmid=None,
                                arxiv_id=None,
                                year=year,
                                journal=item.get("container-title", [""])[0],
                                url=item.get("URL", f"https://doi.org/{item.get('DOI')}"),
                                source=SearchSource.CROSSREF,
                                citations=item.get("is-referenced-by-count", 0)
                            ))
                            
                        except Exception:
                            continue
                            
        except Exception:
            pass
            
        return papers
        
    async def _search_semantic_scholar(self, query: str, max_results: int) -> List[Paper]:
        """搜索Semantic Scholar"""
        papers = []
        
        try:
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                # Semantic Scholar API
                url = f"https://api.semanticscholar.org/graph/v1/paper/search"
                params = {
                    "query": query,
                    "limit": max_results,
                    "fields": "title,authors,abstract,year,venue,citationCount,url,externalIds"
                }
                
                response = await client.get(url, params=params)
                
                if response.status_code == 200:
                    data = response.json()
                    
                    for item in data.get("data", []):
                        try:
                            authors = [a.get("name", "") for a in item.get("authors", [])[:5]]
                            
                            # 提取DOI
                            external_ids = item.get("externalIds", {})
                            doi = external_ids.get("DOI")
                            pmid = external_ids.get("PubMed")
                            arxiv_id = external_ids.get("ArXiv")
                            
                            papers.append(Paper(
                                title=item.get("title", ""),
                                authors=authors,
                                abstract=item.get("abstract", ""),
                                doi=doi,
                                pmid=pmid,
                                arxiv_id=arxiv_id,
                                year=item.get("year", datetime.now().year),
                                journal=item.get("venue", ""),
                                url=item.get("url", ""),
                                source=SearchSource.SEMANTIC_SCHOLAR,
                                citations=item.get("citationCount", 0)
                            ))
                            
                        except Exception:
                            continue
                            
        except Exception:
            pass
            
        return papers
        
    def _deduplicate_papers(self, papers: List[Paper]) -> List[Paper]:
        """去重论文（基于标题相似度和DOI）"""
        seen_dois = set()
        seen_titles = set()
        unique = []
        
        for paper in papers:
            # DOI去重
            if paper.doi and paper.doi in seen_dois:
                continue
                
            # 标题去重（简单处理）
            title_key = re.sub(r'[^a-z0-9]', '', paper.title.lower())[:50]
            if title_key in seen_titles:
                continue
                
            if paper.doi:
                seen_dois.add(paper.doi)
            seen_titles.add(title_key)
            unique.append(paper)
            
        return unique
        
    def _calculate_relevance(self, paper: Paper, query: str) -> float:
        """计算论文相关性分数"""
        score = 0.0
        query_lower = query.lower()
        query_terms = set(query_lower.split())
        
        # 标题匹配（权重最高）
        title_lower = paper.title.lower()
        for term in query_terms:
            if term in title_lower:
                score += 3.0
                
        # 摘要匹配
        abstract_lower = paper.abstract.lower()
        for term in query_terms:
            if term in abstract_lower:
                score += 1.0
                
        # 关键词匹配
        if paper.keywords:
            keywords_lower = [kw.lower() for kw in paper.keywords]
            for term in query_terms:
                if any(term in kw for kw in keywords_lower):
                    score += 2.0
                    
        # 引用数加成
        if paper.citations:
            score += min(paper.citations / 100, 5.0)  # 最多加5分
            
        # 年份加成（越新越好）
        years_old = datetime.now().year - paper.year
        if years_old <= 2:
            score += 2.0
        elif years_old <= 5:
            score += 1.0
            
        return score
        
    def _calculate_statistics(self, papers: List[Paper]) -> Dict[str, Any]:
        """计算统计信息"""
        if not papers:
            return {}
            
        years = [p.year for p in papers]
        citations = [p.citations for p in papers if p.citations is not None]
        
        stats = {
            "year_range": f"{min(years)}-{max(years)}",
            "average_year": sum(years) / len(years),
            "sources": {},
            "top_journals": {},
            "top_authors": {}
        }
        
        # 统计来源
        for paper in papers:
            source = paper.source.value
            stats["sources"][source] = stats["sources"].get(source, 0) + 1
            
            # 统计期刊
            if paper.journal:
                stats["top_journals"][paper.journal] = stats["top_journals"].get(paper.journal, 0) + 1
                
            # 统计作者
            for author in paper.authors:
                stats["top_authors"][author] = stats["top_authors"].get(author, 0) + 1
                
        # 只保留前5个
        stats["top_journals"] = dict(sorted(stats["top_journals"].items(), 
                                          key=lambda x: x[1], reverse=True)[:5])
        stats["top_authors"] = dict(sorted(stats["top_authors"].items(), 
                                         key=lambda x: x[1], reverse=True)[:5])
        
        # 引用统计
        if citations:
            stats["citations"] = {
                "total": sum(citations),
                "average": sum(citations) / len(citations),
                "max": max(citations),
                "min": min(citations)
            }
            
        return stats
        
    def _get_cache_key(self, query: str, sources: List[SearchSource], 
                      year_from: Optional[int], year_to: Optional[int]) -> str:
        """生成缓存键"""
        key_parts = [
            query,
            ",".join([s.value for s in sources]),
            str(year_from or ""),
            str(year_to or "")
        ]
        key_string = "|".join(key_parts)
        return hashlib.md5(key_string.encode()).hexdigest()
