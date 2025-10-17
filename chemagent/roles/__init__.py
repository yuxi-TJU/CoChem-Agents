"""
Chemistry Roles - Markdown-based like SuperClaude_Framework
All roles are defined in markdown files in the /roles directory
"""

from .loader import (
    role_loader,
    get_role,
    list_available_roles
)

__all__ = [
    'role_loader',
    'get_role', 
    'list_available_roles'
]