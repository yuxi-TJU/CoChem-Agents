"""
Base enhancer class for AI coding assistants
"""

from typing import Dict, Any, List, Optional
from abc import ABC, abstractmethod
from pathlib import Path
import json
import yaml


class BaseEnhancer(ABC):
    """Base class for AI assistant enhancers"""
    
    def __init__(self):
        self.commands = {}
        self.sub_agents = {}
        self.tools = {}
        self.prompts = {}
        self.config = self.load_config()
        
    @abstractmethod
    def get_name(self) -> str:
        """Get enhancer name"""
        pass
    
    @abstractmethod
    def get_platform(self) -> str:
        """Get target platform (claude-code, gemini-cli, etc.)"""
        pass
    
    def load_config(self) -> Dict[str, Any]:
        """Load enhancer configuration"""
        config_path = Path.home() / ".chemagent" / f"{self.get_platform()}_config.yaml"
        if config_path.exists():
            with open(config_path, "r") as f:
                return yaml.safe_load(f)
        return self.get_default_config()
    
    @abstractmethod
    def get_default_config(self) -> Dict[str, Any]:
        """Get default configuration"""
        pass
    
    def register_command(self, command):
        """Register a command"""
        self.commands[command.get_name()] = command
        for alias in command.get_aliases():
            self.commands[alias] = command
    
    def register_sub_agent(self, name: str, agent):
        """Register a sub-agent"""
        self.sub_agents[name] = agent
    
    def register_tool(self, name: str, tool):
        """Register a tool"""
        self.tools[name] = tool
    
    def register_prompt(self, name: str, prompt: str):
        """Register a prompt template"""
        self.prompts[name] = prompt
    
    @abstractmethod
    def install(self) -> bool:
        """Install the enhancer"""
        pass
    
    @abstractmethod
    def uninstall(self) -> bool:
        """Uninstall the enhancer"""
        pass
    
    def get_command(self, name: str):
        """Get a registered command"""
        return self.commands.get(name)
    
    def list_commands(self) -> List[str]:
        """List all registered commands"""
        return list(self.commands.keys())
    
    def list_sub_agents(self) -> List[str]:
        """List all registered sub-agents"""
        return list(self.sub_agents.keys())
    
    def list_tools(self) -> List[str]:
        """List all registered tools"""
        return list(self.tools.keys())
    
    async def execute_command(self, command_name: str, *args, **kwargs) -> Any:
        """Execute a command"""
        command = self.get_command(command_name)
        if command:
            return await command.execute(*args, **kwargs)
        raise ValueError(f"Command not found: {command_name}")
    
    def get_help(self) -> str:
        """Get help text for the enhancer"""
        help_text = f"""
# {self.get_name()} - Chemistry Enhancement for {self.get_platform()}

## Available Commands:
"""
        for name, cmd in self.commands.items():
            if not any(name == alias for alias in cmd.get_aliases()):
                help_text += f"- **{name}**: {cmd.get_description()}\n"
        
        help_text += "\n## Sub-Agents:\n"
        for name, agent in self.sub_agents.items():
            help_text += f"- **@{name}**: {agent.get_description()}\n"
        
        help_text += "\n## Tools:\n"
        for name in self.tools.keys():
            help_text += f"- **{name}**\n"
        
        return help_text
