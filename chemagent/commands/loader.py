"""
Markdown-based Command Loader
Similar to SuperClaude_Framework's approach
"""

import os
from pathlib import Path
from typing import Dict, Any, List, Optional
import yaml
import re


class CommandLoader:
    """
    Load commands from Markdown files
    Following SuperClaude_Framework pattern
    """
    
    def __init__(self):
        # Command directories (in priority order)
        self.command_dirs = [
            Path(".claude/commands"),  # Project-specific commands
            Path.home() / ".claude/commands",  # Global user commands
            Path(__file__).parent.parent.parent / "commands",  # Package commands
        ]
        
        self.commands = {}
        self.load_all_commands()
    
    def load_all_commands(self):
        """Load all available commands from markdown files"""
        for cmd_dir in self.command_dirs:
            if cmd_dir.exists() and cmd_dir.is_dir():
                for md_file in cmd_dir.glob("*.md"):
                    self.load_command(md_file)
    
    def load_command(self, filepath: Path):
        """Load a single command from markdown file"""
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Extract command name from filename
            cmd_name = filepath.stem
            
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
            
            # Store command
            self.commands[cmd_name] = {
                'name': cmd_name,
                'description': metadata.get('description', f'{cmd_name} command'),
                'content': content,
                'metadata': metadata,
                'filepath': str(filepath),
                'tools': metadata.get('tools', []),
                'models': metadata.get('models', []),
                'parameters': metadata.get('parameters', {}),
            }
            
        except Exception as e:
            print(f"Error loading command {filepath}: {e}")
    
    def get_command(self, name: str) -> Optional[Dict[str, Any]]:
        """Get a specific command by name"""
        return self.commands.get(name)
    
    def list_commands(self) -> List[str]:
        """List all available command names"""
        return list(self.commands.keys())
    
    def get_command_help(self, name: str) -> str:
        """Get help text for a command"""
        cmd = self.get_command(name)
        if not cmd:
            return f"Command '{name}' not found"
        
        help_text = f"# {cmd['name']}\n\n"
        help_text += f"{cmd['description']}\n\n"
        
        # Add parameter info if available
        if cmd['parameters']:
            help_text += "## Parameters\n"
            for param, info in cmd['parameters'].items():
                required = "required" if info.get('required') else "optional"
                help_text += f"- **{param}** ({required}): {info.get('description', '')}\n"
            help_text += "\n"
        
        # Add tools info if specified
        if cmd['tools']:
            help_text += f"## Required Tools\n"
            help_text += f"{', '.join(cmd['tools'])}\n\n"
        
        return help_text
    
    def execute_command(self, name: str, args: List[str] = None, kwargs: Dict[str, Any] = None) -> str:
        """
        Get command prompt for execution
        This returns the markdown content that should be sent to the AI
        """
        cmd = self.get_command(name)
        if not cmd:
            return f"Error: Command '{name}' not found. Available commands: {', '.join(self.list_commands())}"
        
        # Get the command content
        prompt = cmd['content']
        
        # Replace parameters if any
        if args:
            # Simple positional parameter replacement
            for i, arg in enumerate(args):
                prompt = prompt.replace(f"${i+1}", str(arg))
                prompt = prompt.replace(f"${{arg{i+1}}}", str(arg))
        
        if kwargs:
            # Named parameter replacement
            for key, value in kwargs.items():
                prompt = prompt.replace(f"${{{key}}}", str(value))
                prompt = prompt.replace(f"${key}", str(value))
        
        # Add context about available tools
        if cmd['tools']:
            prompt += f"\n\n**Note**: This command uses the following tools: {', '.join(cmd['tools'])}"
        
        return prompt
    
    def create_command(self, name: str, description: str, content: str, 
                      tools: List[str] = None, save_location: str = "project") -> bool:
        """
        Create a new command and save it as markdown
        
        Args:
            name: Command name (will be the filename)
            description: Command description
            content: Command prompt/instructions
            tools: List of required tools
            save_location: "project" or "global"
        """
        # Determine save directory
        if save_location == "project":
            save_dir = Path(".claude/commands")
        else:
            save_dir = Path.home() / ".claude/commands"
        
        # Create directory if needed
        save_dir.mkdir(parents=True, exist_ok=True)
        
        # Build markdown content
        md_content = "---\n"
        md_content += f"description: {description}\n"
        if tools:
            md_content += f"tools: {tools}\n"
        md_content += "---\n\n"
        md_content += content
        
        # Save file
        filepath = save_dir / f"{name}.md"
        try:
            with open(filepath, 'w', encoding='utf-8') as f:
                f.write(md_content)
            
            # Reload commands
            self.load_command(filepath)
            return True
            
        except Exception as e:
            print(f"Error saving command: {e}")
            return False


# Global command loader instance
command_loader = CommandLoader()


def register_command(name: str, description: str, content: str, **kwargs):
    """
    Register a new command dynamically
    Similar to SuperClaude's approach
    """
    return command_loader.create_command(name, description, content, **kwargs)


def get_command_prompt(name: str, *args, **kwargs) -> str:
    """
    Get the prompt for a command with parameters filled in
    This is what gets sent to the AI
    """
    return command_loader.execute_command(name, args, kwargs)


def list_available_commands() -> List[Dict[str, str]]:
    """List all available commands with descriptions"""
    results = []
    for cmd_name in command_loader.list_commands():
        cmd = command_loader.get_command(cmd_name)
        results.append({
            'name': cmd_name,
            'description': cmd['description']
        })
    return results
