"""
ChemAgent Commands - Markdown-based like SuperClaude_Framework
"""

from .loader import (
    command_loader,
    get_command_prompt,
    list_available_commands,
    register_command
)

# Export the main interface
__all__ = [
    'command_loader',
    'get_command_prompt',
    'list_available_commands',
    'register_command',
]

# Initialize and load all commands on import
# This will load from:
# 1. .claude/commands/ (project-specific)
# 2. ~/.claude/commands/ (user global)
# 3. Package commands directory
command_loader.load_all_commands()