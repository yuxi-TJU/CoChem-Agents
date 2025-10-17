#!/usr/bin/env python3
"""
Enhanced Base MCP Server Implementation (v2)
增强版 MCP 服务器基类 - 更符合官方规范
"""

import asyncio
import json
import logging
import sys
from typing import Any, Dict, List, Optional, Union
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class RiskLevel(Enum):
    """Risk levels for tool annotations"""
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"
    CRITICAL = "critical"


@dataclass
class ToolAnnotation:
    """Tool annotation for risk and permission management"""
    risk_level: RiskLevel = RiskLevel.LOW
    requires_confirmation: bool = False
    modifies_data: bool = False
    accesses_network: bool = False
    accesses_filesystem: bool = False
    description: str = ""


@dataclass
class MCPRequest:
    """MCP Request structure (JSON-RPC 2.0)"""
    jsonrpc: str = "2.0"
    method: str = ""
    params: Union[Dict[str, Any], List, None] = None
    id: Optional[Union[str, int]] = None


@dataclass
class MCPResponse:
    """MCP Response structure (JSON-RPC 2.0)"""
    jsonrpc: str = "2.0"
    result: Optional[Any] = None
    error: Optional[Dict[str, Any]] = None
    id: Optional[Union[str, int]] = None


@dataclass
class MCPError:
    """Standard JSON-RPC 2.0 error"""
    code: int
    message: str
    data: Optional[Any] = None
    
    # Standard error codes
    PARSE_ERROR = -32700
    INVALID_REQUEST = -32600
    METHOD_NOT_FOUND = -32601
    INVALID_PARAMS = -32602
    INTERNAL_ERROR = -32603
    
    def to_dict(self) -> Dict[str, Any]:
        result = {"code": self.code, "message": self.message}
        if self.data is not None:
            result["data"] = self.data
        return result


class BaseMCPServerV2(ABC):
    """
    Enhanced base class for MCP servers
    Implements Model Context Protocol with better compliance
    """
    
    # Protocol version
    PROTOCOL_VERSION = "1.0.0"
    
    def __init__(self, name: str, version: str = "1.0.0"):
        self.name = name
        self.version = version
        self.capabilities = {}
        self.tools = {}
        self.tool_annotations = {}
        self.client_info = {}
        self.session_id = None
        self._setup_capabilities()
        self._register_tools()
        
    @abstractmethod
    def _setup_capabilities(self):
        """Setup server capabilities - must be implemented by subclasses"""
        pass
        
    @abstractmethod
    def _register_tools(self):
        """Register available tools - must be implemented by subclasses"""
        pass
    
    def get_server_info(self) -> Dict[str, Any]:
        """Get server information for initialize response"""
        return {
            "protocolVersion": self.PROTOCOL_VERSION,
            "serverInfo": {
                "name": self.name,
                "version": self.version
            },
            "capabilities": self.capabilities,
            "sessionId": self.session_id or f"{self.name}-{datetime.now().isoformat()}"
        }
    
    async def handle_initialize(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Handle initialize request with version negotiation"""
        # Store client information
        self.client_info = params.get("clientInfo", {})
        client_version = params.get("protocolVersion", "1.0.0")
        
        # Version negotiation (simplified)
        if client_version != self.PROTOCOL_VERSION:
            logger.warning(f"Client version {client_version} != Server version {self.PROTOCOL_VERSION}")
        
        # Generate session ID
        self.session_id = f"{self.name}-{datetime.now().isoformat()}"
        
        return self.get_server_info()
    
    async def handle_list_tools(self, params: Dict[str, Any]) -> List[Dict[str, Any]]:
        """List available tools with annotations"""
        tools = []
        
        for tool_name, tool_func in self.tools.items():
            tool_info = {
                "name": tool_name,
                "description": tool_func.__doc__ or "No description",
                "inputSchema": getattr(tool_func, "_mcp_params", {})
            }
            
            # Add annotations if available
            if tool_name in self.tool_annotations:
                annotation = self.tool_annotations[tool_name]
                tool_info["annotations"] = {
                    "riskLevel": annotation.risk_level.value,
                    "requiresConfirmation": annotation.requires_confirmation,
                    "modifiesData": annotation.modifies_data,
                    "accessesNetwork": annotation.accesses_network,
                    "accessesFilesystem": annotation.accesses_filesystem
                }
            
            tools.append(tool_info)
        
        return tools
    
    async def handle_tool_call(self, tool_name: str, arguments: Dict[str, Any]) -> Any:
        """Execute a tool with the given arguments"""
        if tool_name not in self.tools:
            raise ValueError(f"Tool not found: {tool_name}")
        
        tool_func = self.tools[tool_name]
        
        # Check annotations for confirmation requirement
        if tool_name in self.tool_annotations:
            annotation = self.tool_annotations[tool_name]
            if annotation.requires_confirmation:
                logger.warning(f"Tool {tool_name} requires confirmation (risk: {annotation.risk_level.value})")
        
        # Execute tool
        if asyncio.iscoroutinefunction(tool_func):
            return await tool_func(**arguments)
        else:
            return tool_func(**arguments)
    
    async def handle_completions(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Handle parameter completions request (for developer experience)"""
        # This is a placeholder - implement based on your tools
        return {
            "completions": [],
            "isIncomplete": False
        }
    
    async def process_request(self, request: MCPRequest) -> MCPResponse:
        """Process a single MCP request"""
        try:
            # Route to appropriate handler
            if request.method == "initialize":
                result = await self.handle_initialize(request.params or {})
            elif request.method == "tools/list":
                result = await self.handle_list_tools(request.params or {})
            elif request.method == "tools/call":
                tool_name = request.params.get("name")
                arguments = request.params.get("arguments", {})
                result = await self.handle_tool_call(tool_name, arguments)
            elif request.method == "completions":
                result = await self.handle_completions(request.params or {})
            else:
                # Try custom method handler
                result = await self.handle_custom_method(request.method, request.params)
            
            return MCPResponse(result=result, id=request.id)
            
        except ValueError as e:
            error = MCPError(
                code=MCPError.INVALID_PARAMS,
                message=str(e)
            )
            return MCPResponse(error=error.to_dict(), id=request.id)
        except NotImplementedError:
            error = MCPError(
                code=MCPError.METHOD_NOT_FOUND,
                message=f"Method not found: {request.method}"
            )
            return MCPResponse(error=error.to_dict(), id=request.id)
        except Exception as e:
            logger.error(f"Error processing request: {e}", exc_info=True)
            error = MCPError(
                code=MCPError.INTERNAL_ERROR,
                message="Internal server error",
                data=str(e)
            )
            return MCPResponse(error=error.to_dict(), id=request.id)
    
    async def handle_custom_method(self, method: str, params: Any) -> Any:
        """Handle custom methods - override in subclasses"""
        raise NotImplementedError(f"Method not implemented: {method}")
    
    async def process_message(self, message: str) -> str:
        """Process incoming JSON-RPC message (supports batching)"""
        try:
            data = json.loads(message)
            
            # Check if it's a batch request
            if isinstance(data, list):
                # Batch request
                responses = []
                for request_data in data:
                    request = MCPRequest(
                        jsonrpc=request_data.get("jsonrpc", "2.0"),
                        method=request_data.get("method"),
                        params=request_data.get("params"),
                        id=request_data.get("id")
                    )
                    response = await self.process_request(request)
                    if request.id is not None:  # Only include responses for requests with IDs
                        responses.append(response.__dict__)
                
                return json.dumps(responses)
            else:
                # Single request
                request = MCPRequest(
                    jsonrpc=data.get("jsonrpc", "2.0"),
                    method=data.get("method"),
                    params=data.get("params"),
                    id=data.get("id")
                )
                response = await self.process_request(request)
                return json.dumps(response.__dict__)
                
        except json.JSONDecodeError as e:
            error = MCPError(
                code=MCPError.PARSE_ERROR,
                message="Parse error",
                data=str(e)
            )
            response = MCPResponse(error=error.to_dict())
            return json.dumps(response.__dict__)
    
    def register_tool(self, name: str, func: callable, annotation: Optional[ToolAnnotation] = None):
        """Register a tool with optional annotations"""
        self.tools[name] = func
        if annotation:
            self.tool_annotations[name] = annotation
    
    async def send_progress(self, progress: Dict[str, Any]):
        """Send progress notification (for long-running operations)"""
        # This would send a notification in a real implementation
        logger.info(f"Progress: {progress}")
    
    async def run_stdio(self):
        """Run server using stdin/stdout (standard MCP mode)"""
        logger.info(f"Starting {self.name} MCP Server v{self.version} (stdio mode)")
        
        reader = asyncio.StreamReader()
        protocol = asyncio.StreamReaderProtocol(reader)
        await asyncio.get_event_loop().connect_read_pipe(
            lambda: protocol, sys.stdin
        )
        
        while True:
            try:
                # Read line from stdin
                line = await reader.readline()
                if not line:
                    break
                    
                message = line.decode().strip()
                if not message:
                    continue
                    
                # Process message
                response = await self.process_message(message)
                
                # Write response to stdout
                print(response)
                sys.stdout.flush()
                
            except KeyboardInterrupt:
                logger.info("Server shutting down...")
                break
            except Exception as e:
                logger.error(f"Error in main loop: {e}")


def mcp_tool(
    params: Dict[str, Any] = None,
    annotation: Optional[ToolAnnotation] = None
):
    """
    Decorator to mark methods as MCP tools with annotations
    
    Example:
        @mcp_tool(
            params={"text": {"type": "string", "required": True}},
            annotation=ToolAnnotation(
                risk_level=RiskLevel.LOW,
                modifies_data=False
            )
        )
        def process_text(self, text: str) -> str:
            return text.upper()
    """
    def decorator(func):
        func._mcp_params = params or {}
        func._mcp_annotation = annotation or ToolAnnotation()
        return func
    return decorator
