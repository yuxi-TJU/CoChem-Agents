"""
Literature search tools for ChemAgent
文献搜索工具实现
"""

from typing import Dict, List, Any, Optional
import logging
import asyncio
import requests
from datetime import datetime

from ..core.registry import BaseTool, ToolMetadata, ToolCategory


logger = logging.getLogger(__name__)


class LiteratureSearchTool(BaseTool):
    """Tool for searching scientific literature"""
    
    def get_metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="literature_search",
            description="Search scientific literature from multiple databases",
            category=ToolCategory.LITERATURE,
            version="1.0.0",
            author="ChemAgent",
            dependencies=["requests", "biopython"],
            parameters={
                "required": ["query"],
                "properties": {
                    "query": {"type": "string", "description": "Search query"},
                    "database": {
                        "type": "string",
                        "enum": ["pubmed", "chemrxiv", "scholar", "crossref", "all"],
                        "description": "Database to search",
                        "default": "pubmed"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of results",
                        "default": 10
                    },
                    "year_from": {
                        "type": "integer",
                        "description": "Start year for filtering"
                    },
                    "year_to": {
                        "type": "integer",
                        "description": "End year for filtering"
                    },
                    "sort_by": {
                        "type": "string",
                        "enum": ["relevance", "date", "citations"],
                        "default": "relevance"
                    }
                }
            },
            examples=[
                {"query": "PROTAC drug design", "database": "pubmed", "limit": 5},
                {"query": "machine learning molecular property", "database": "all"},
                {"query": "organic synthesis green chemistry", "year_from": 2020}
            ]
        )
    
    async def execute(self,
                      query: str,
                      database: str = "pubmed",
                      limit: int = 10,
                      year_from: Optional[int] = None,
                      year_to: Optional[int] = None,
                      sort_by: str = "relevance",
                      **kwargs) -> Dict[str, Any]:
        """Execute literature search"""
        try:
            results = {
                "query": query,
                "database": database,
                "total_results": 0,
                "papers": []
            }
            
            if database == "pubmed" or database == "all":
                pubmed_results = await self._search_pubmed(
                    query, limit, year_from, year_to
                )
                results["papers"].extend(pubmed_results)
            
            if database == "chemrxiv" or database == "all":
                chemrxiv_results = await self._search_chemrxiv(
                    query, limit
                )
                results["papers"].extend(chemrxiv_results)
            
            if database == "crossref" or database == "all":
                crossref_results = await self._search_crossref(
                    query, limit, year_from, year_to
                )
                results["papers"].extend(crossref_results)
            
            if database == "scholar" or database == "all":
                # Google Scholar 需要特殊处理（可能需要代理）
                scholar_results = await self._search_scholar(
                    query, limit
                )
                results["papers"].extend(scholar_results)
            
            # 排序结果
            if sort_by == "date":
                results["papers"].sort(
                    key=lambda x: x.get("year", 0),
                    reverse=True
                )
            elif sort_by == "citations":
                results["papers"].sort(
                    key=lambda x: x.get("citations", 0),
                    reverse=True
                )
            
            # 限制结果数量
            if database == "all" and len(results["papers"]) > limit:
                results["papers"] = results["papers"][:limit]
            
            results["total_results"] = len(results["papers"])
            
            return results
            
        except Exception as e:
            logger.error(f"Literature search failed: {e}")
            return {"error": str(e)}
    
    async def _search_pubmed(self,
                            query: str,
                            limit: int,
                            year_from: Optional[int] = None,
                            year_to: Optional[int] = None) -> List[Dict[str, Any]]:
        """Search PubMed database"""
        try:
            from Bio import Entrez
            Entrez.email = "chemagent@example.com"  # Required by NCBI
            
            # Build query with date filters
            search_query = query
            if year_from or year_to:
                year_from = year_from or 1900
                year_to = year_to or datetime.now().year
                search_query += f" AND {year_from}:{year_to}[dp]"
            
            # Search for IDs
            handle = Entrez.esearch(
                db="pubmed",
                term=search_query,
                retmax=limit,
                sort="relevance"
            )
            record = Entrez.read(handle)
            handle.close()
            
            ids = record.get("IdList", [])
            if not ids:
                return []
            
            # Fetch details
            handle = Entrez.efetch(
                db="pubmed",
                id=ids,
                rettype="medline",
                retmode="xml"
            )
            records = Entrez.read(handle)
            handle.close()
            
            papers = []
            for article in records.get("PubmedArticle", []):
                medline = article.get("MedlineCitation", {})
                article_data = medline.get("Article", {})
                
                # Extract authors
                authors = []
                author_list = article_data.get("AuthorList", [])
                for author in author_list[:3]:  # First 3 authors
                    last = author.get("LastName", "")
                    first = author.get("ForeName", "")
                    if last:
                        authors.append(f"{last} {first[0] if first else ''}")
                
                # Extract publication year
                pub_date = medline.get("DateCompleted", {})
                year = pub_date.get("Year", "")
                
                paper = {
                    "title": article_data.get("ArticleTitle", ""),
                    "authors": authors,
                    "journal": article_data.get("Journal", {}).get("Title", ""),
                    "year": int(year) if year else None,
                    "pmid": medline.get("PMID", ""),
                    "abstract": article_data.get("Abstract", {}).get("AbstractText", [""])[0],
                    "doi": self._extract_doi(article),
                    "database": "PubMed"
                }
                papers.append(paper)
            
            return papers
            
        except Exception as e:
            logger.error(f"PubMed search failed: {e}")
            return []
    
    async def _search_chemrxiv(self, query: str, limit: int) -> List[Dict[str, Any]]:
        """Search ChemRxiv preprints"""
        try:
            # ChemRxiv API endpoint
            api_url = "https://chemrxiv.org/engage/chemrxiv/public-api/v1/items"
            
            params = {
                "term": query,
                "limit": limit,
                "skip": 0
            }
            
            response = requests.get(api_url, params=params, timeout=10)
            if response.status_code != 200:
                return []
            
            data = response.json()
            papers = []
            
            for item in data.get("itemHits", []):
                paper = {
                    "title": item.get("item", {}).get("title", ""),
                    "authors": self._extract_chemrxiv_authors(item),
                    "year": self._extract_year(item.get("item", {}).get("publishedDate")),
                    "abstract": item.get("item", {}).get("abstract", ""),
                    "doi": item.get("item", {}).get("doi", ""),
                    "database": "ChemRxiv",
                    "preprint": True
                }
                papers.append(paper)
            
            return papers
            
        except Exception as e:
            logger.error(f"ChemRxiv search failed: {e}")
            return []
    
    async def _search_crossref(self,
                               query: str,
                               limit: int,
                               year_from: Optional[int] = None,
                               year_to: Optional[int] = None) -> List[Dict[str, Any]]:
        """Search CrossRef database"""
        try:
            api_url = "https://api.crossref.org/works"
            
            params = {
                "query": query,
                "rows": limit
            }
            
            # Add date filters
            if year_from:
                params["filter"] = f"from-pub-date:{year_from}"
            if year_to:
                if "filter" in params:
                    params["filter"] += f",until-pub-date:{year_to}"
                else:
                    params["filter"] = f"until-pub-date:{year_to}"
            
            response = requests.get(api_url, params=params, timeout=10)
            if response.status_code != 200:
                return []
            
            data = response.json()
            papers = []
            
            for item in data.get("message", {}).get("items", []):
                # Extract authors
                authors = []
                for author in item.get("author", [])[:3]:
                    given = author.get("given", "")
                    family = author.get("family", "")
                    if family:
                        authors.append(f"{family} {given[0] if given else ''}")
                
                paper = {
                    "title": item.get("title", [""])[0] if item.get("title") else "",
                    "authors": authors,
                    "journal": item.get("container-title", [""])[0] if item.get("container-title") else "",
                    "year": item.get("published-print", {}).get("date-parts", [[None]])[0][0],
                    "doi": item.get("DOI", ""),
                    "citations": item.get("is-referenced-by-count", 0),
                    "database": "CrossRef"
                }
                papers.append(paper)
            
            return papers
            
        except Exception as e:
            logger.error(f"CrossRef search failed: {e}")
            return []
    
    async def _search_scholar(self, query: str, limit: int) -> List[Dict[str, Any]]:
        """Search Google Scholar (limited functionality without API)"""
        try:
            # Note: Google Scholar doesn't have official API
            # This is a simplified implementation
            # For production, consider using scholarly library or scraping carefully
            
            papers = []
            
            # Placeholder for Google Scholar integration
            # In production, you might use:
            # - scholarly library
            # - SerpAPI (paid)
            # - Custom scraping (be careful with rate limits)
            
            logger.info("Google Scholar search not fully implemented")
            
            return papers
            
        except Exception as e:
            logger.error(f"Google Scholar search failed: {e}")
            return []
    
    def _extract_doi(self, pubmed_article: Dict) -> str:
        """Extract DOI from PubMed article"""
        article_ids = pubmed_article.get("PubmedData", {}).get("ArticleIdList", [])
        for aid in article_ids:
            if aid.attributes.get("IdType") == "doi":
                return str(aid)
        return ""
    
    def _extract_chemrxiv_authors(self, item: Dict) -> List[str]:
        """Extract authors from ChemRxiv item"""
        authors = []
        author_list = item.get("item", {}).get("authors", [])
        for author in author_list[:3]:
            name = author.get("firstName", "") + " " + author.get("lastName", "")
            authors.append(name.strip())
        return authors
    
    def _extract_year(self, date_str: str) -> Optional[int]:
        """Extract year from date string"""
        if not date_str:
            return None
        try:
            return int(date_str[:4])
        except:
            return None


class CitationAnalyzer(BaseTool):
    """Tool for analyzing citations and finding related papers"""
    
    def get_metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="citation_analyzer",
            description="Analyze citations and find related papers",
            category=ToolCategory.LITERATURE,
            version="1.0.0",
            author="ChemAgent",
            dependencies=["requests"],
            parameters={
                "required": ["doi_or_pmid"],
                "properties": {
                    "doi_or_pmid": {
                        "type": "string",
                        "description": "DOI or PubMed ID of the paper"
                    },
                    "analysis_type": {
                        "type": "string",
                        "enum": ["citations", "references", "related", "metrics"],
                        "default": "related"
                    }
                }
            },
            examples=[
                {"doi_or_pmid": "10.1038/nature12373", "analysis_type": "citations"},
                {"doi_or_pmid": "35524985", "analysis_type": "related"}
            ]
        )
    
    async def execute(self,
                      doi_or_pmid: str,
                      analysis_type: str = "related",
                      **kwargs) -> Dict[str, Any]:
        """Analyze citations for a paper"""
        try:
            result = {
                "input": doi_or_pmid,
                "analysis_type": analysis_type,
                "data": {}
            }
            
            if analysis_type == "citations":
                # Get papers that cite this paper
                result["data"]["citing_papers"] = await self._get_citations(doi_or_pmid)
                
            elif analysis_type == "references":
                # Get papers referenced by this paper
                result["data"]["references"] = await self._get_references(doi_or_pmid)
                
            elif analysis_type == "related":
                # Get related papers
                result["data"]["related_papers"] = await self._get_related(doi_or_pmid)
                
            elif analysis_type == "metrics":
                # Get citation metrics
                result["data"]["metrics"] = await self._get_metrics(doi_or_pmid)
            
            return result
            
        except Exception as e:
            logger.error(f"Citation analysis failed: {e}")
            return {"error": str(e)}
    
    async def _get_citations(self, doi_or_pmid: str) -> List[Dict[str, Any]]:
        """Get papers that cite the given paper"""
        # Implementation would use CrossRef or other citation databases
        return []
    
    async def _get_references(self, doi_or_pmid: str) -> List[Dict[str, Any]]:
        """Get references of the given paper"""
        # Implementation would parse paper references
        return []
    
    async def _get_related(self, doi_or_pmid: str) -> List[Dict[str, Any]]:
        """Get related papers"""
        # Implementation would use similarity algorithms
        return []
    
    async def _get_metrics(self, doi_or_pmid: str) -> Dict[str, Any]:
        """Get citation metrics"""
        # Implementation would calculate h-index, impact factor, etc.
        return {
            "citation_count": 0,
            "h_index": 0,
            "altmetric_score": 0
        }
