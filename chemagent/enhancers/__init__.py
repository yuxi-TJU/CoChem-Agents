"""
Enhancers for different AI coding assistants
Provides chemistry capabilities to Claude Code, Gemini CLI, etc.
"""

from .claude_enhancer import ClaudeCodeEnhancer
from .gemini_enhancer import GeminiCLIEnhancer
from .base import BaseEnhancer

__all__ = [
    "ClaudeCodeEnhancer",
    "GeminiCLIEnhancer",
    "BaseEnhancer",
]
