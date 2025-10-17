"""ChemAgent - Chemistry Enhancement Package for AI Assistants"""

__version__ = "0.1.0"

from .core.config import ChemAgentConfig
from .enhancers import ClaudeCodeEnhancer, GeminiCLIEnhancer
from .commands import list_available_commands, get_command_prompt
from .roles import list_available_roles, get_role

__all__ = [
    "ChemAgentConfig",
    "ClaudeCodeEnhancer",
    "GeminiCLIEnhancer",
    "list_available_commands",
    "get_command_prompt",
    "list_available_roles",
    "get_role",
]
