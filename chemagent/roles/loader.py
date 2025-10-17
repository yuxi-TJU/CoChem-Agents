"""
Markdown-based Role Loader
Similar to SuperClaude_Framework's persona system
"""

from pathlib import Path
from typing import Dict, Any, List, Optional
import yaml


class RoleLoader:
    """
    Load roles/personas from Markdown files
    Following SuperClaude_Framework pattern
    """
    
    def __init__(self):
        # Role directories (in priority order)
        self.role_dirs = [
            Path(".claude/roles"),  # Project-specific roles
            Path.home() / ".claude/roles",  # Global user roles
            Path(__file__).parent.parent.parent / "roles",  # Package roles
        ]
        
        self.roles = {}
        self.load_all_roles()
    
    def load_all_roles(self):
        """Load all available roles from markdown files"""
        for role_dir in self.role_dirs:
            if role_dir.exists() and role_dir.is_dir():
                for md_file in role_dir.glob("*.md"):
                    self.load_role(md_file)
    
    def load_role(self, filepath: Path):
        """Load a single role from markdown file"""
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Extract role name from filename
            role_name = filepath.stem
            
            # Parse YAML frontmatter if present
            metadata = {}
            if content.startswith('---'):
                parts = content.split('---', 2)
                if len(parts) >= 3:
                    try:
                        metadata = yaml.safe_load(parts[1])
                    except:
                        pass
                    content = parts[2].strip()
            
            # Store role
            self.roles[role_name] = {
                'name': metadata.get('name', role_name),
                'description': metadata.get('description', f'{role_name} role'),
                'content': content,
                'metadata': metadata,
                'filepath': str(filepath),
                'specialties': metadata.get('specialties', []),
                'tools': metadata.get('tools', []),
            }
            
        except Exception as e:
            print(f"Error loading role {filepath}: {e}")
    
    def get_role(self, name: str) -> Optional[Dict[str, Any]]:
        """Get a specific role by name"""
        # Try exact match first
        if name in self.roles:
            return self.roles[name]
        
        # Try with @ prefix removed
        if name.startswith('@'):
            name_without_at = name[1:]
            if name_without_at in self.roles:
                return self.roles[name_without_at]
        
        # Try fuzzy match
        for role_name in self.roles:
            if name.lower() in role_name.lower():
                return self.roles[role_name]
        
        return None
    
    def list_roles(self) -> List[str]:
        """List all available role names"""
        return list(self.roles.keys())
    
    def get_role_prompt(self, name: str) -> str:
        """Get the role prompt/persona definition"""
        role = self.get_role(name)
        if not role:
            return f"Role '{name}' not found. Available roles: {', '.join(self.list_roles())}"
        
        return role['content']


# Global role loader instance
role_loader = RoleLoader()


def get_role(name: str) -> Optional[Dict[str, Any]]:
    """Get role information"""
    return role_loader.get_role(name)


def list_available_roles() -> List[Dict[str, str]]:
    """List all available roles with descriptions"""
    results = []
    for role_name in role_loader.list_roles():
        role = role_loader.get_role(role_name)
        results.append({
            'name': f"@{role_name}",
            'description': role['description']
        })
    return results
