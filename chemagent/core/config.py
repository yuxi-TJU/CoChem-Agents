"""Configuration for ChemAgent Enhancement Package"""

from typing import Dict, List, Optional, Any
from pathlib import Path
from pydantic import BaseModel, Field
from pydantic_settings import BaseSettings
import yaml
import json


class ToolConfig(BaseModel):
    """Tool configuration"""
    name: str = Field(description="Tool name")
    enabled: bool = Field(default=True, description="Whether tool is enabled")
    config: Dict[str, Any] = Field(default_factory=dict, description="Tool-specific configuration")
    

class ChemAgentConfig(BaseSettings):
    """Configuration for ChemAgent Enhancement Package"""
    
    # MCP Server settings
    mcp_host: str = Field(default="localhost", description="MCP server host")
    mcp_port: int = Field(default=8765, description="MCP server port")
    
    # Tool configurations
    tools: List[ToolConfig] = Field(default_factory=list)
    
    # General settings
    workspace_dir: Path = Field(default=Path.home() / ".chemagent", description="Workspace directory")
    cache_dir: Path = Field(default=Path.home() / ".chemagent" / "cache", description="Cache directory")
    log_level: str = Field(default="INFO", description="Logging level")
    
    class Config:
        env_prefix = "CHEMAGENT_"
        env_file = ".env"
        env_file_encoding = "utf-8"
        
    @classmethod
    def from_yaml(cls, path: Path) -> "ChemAgentConfig":
        """Load configuration from YAML file"""
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        return cls(**data)
    
    def ensure_directories(self) -> None:
        """Ensure all required directories exist"""
        self.workspace_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
