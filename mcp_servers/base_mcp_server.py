#!/usr/bin/env python3
"""
Base MCP Server Implementation
基础 MCP 服务器实现，其他化学 MCP 服务器可以继承此类
"""

import asyncio
import json
import logging
from typing import Any, Dict, List, Optional
from abc import ABC, abstractmethod
import sys
from dataclasses import dataclass

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class MCPRequest:
    """MCP Request structure"""
    method: str
    params: Dict[str, Any]
    id: Optional[str] = None


@dataclass
class MCPResponse:
    """MCP Response structure"""
    result: Optional[Any] = None
    error: Optional[Dict[str, Any]] = None
    id: Optional[str] = None


class BaseMCPServer(ABC):
    """
    Base class for MCP servers
    Implements the Model Context Protocol
    """
    
    def __init__(self, name: str, version: str = "1.0.0"):
        self.name = name
        self.version = version
        self.capabilities = {}
        self._setup_capabilities()
        
    @abstractmethod
    def _setup_capabilities(self):
        """Setup server capabilities - must be implemented by subclasses"""
        pass
        
    @abstractmethod
    async def handle_request(self, request: MCPRequest) -> MCPResponse:
        """Handle incoming MCP requests - must be implemented by subclasses"""
        pass
    
    def get_server_info(self) -> Dict[str, Any]:
        """Get server information"""
        return {
            "name": self.name,
            "version": self.version,
            "protocol": "mcp/1.0",
            "capabilities": self.capabilities
        }
    
    async def list_tools(self) -> List[Dict[str, Any]]:
        """List available tools/methods"""
        tools = []
        for method_name in dir(self):
            if method_name.startswith("tool_"):
                method = getattr(self, method_name)
                if callable(method) and hasattr(method, "__doc__"):
                    tools.append({
                        "name": method_name[5:],  # Remove "tool_" prefix
                        "description": method.__doc__ or "No description",
                        "parameters": getattr(method, "_mcp_params", {})
                    })
        return tools
    
    async def process_message(self, message: str) -> str:
        """Process incoming JSON message"""
        try:
            data = json.loads(message)
            request = MCPRequest(
                method=data.get("method"),
                params=data.get("params", {}),
                id=data.get("id")
            )
            
            # Handle built-in methods
            if request.method == "initialize":
                response = MCPResponse(
                    result=self.get_server_info(),
                    id=request.id
                )
            elif request.method == "list_tools":
                response = MCPResponse(
                    result=await self.list_tools(),
                    id=request.id
                )
            else:
                # Delegate to subclass
                response = await self.handle_request(request)
                
            return json.dumps({
                "jsonrpc": "2.0",
                "result": response.result,
                "error": response.error,
                "id": response.id
            })
            
        except json.JSONDecodeError as e:
            return json.dumps({
                "jsonrpc": "2.0",
                "error": {
                    "code": -32700,
                    "message": "Parse error",
                    "data": str(e)
                },
                "id": None
            })
        except Exception as e:
            logger.error(f"Error processing message: {e}")
            return json.dumps({
                "jsonrpc": "2.0",
                "error": {
                    "code": -32603,
                    "message": "Internal error",
                    "data": str(e)
                },
                "id": request.id if 'request' in locals() else None
            })
    
    async def run_stdio(self):
        """Run server using stdin/stdout (standard MCP mode)"""
        logger.info(f"Starting {self.name} MCP Server v{self.version}")
        
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
    
    async def run_tcp(self, host: str = "localhost", port: int = 8765):
        """Run server as TCP server (for testing)"""
        logger.info(f"Starting {self.name} MCP Server v{self.version} on {host}:{port}")
        
        async def handle_client(reader, writer):
            addr = writer.get_extra_info('peername')
            logger.info(f"Client connected from {addr}")
            
            try:
                while True:
                    data = await reader.readline()
                    if not data:
                        break
                    
                    message = data.decode().strip()
                    response = await self.process_message(message)
                    
                    writer.write((response + "\n").encode())
                    await writer.drain()
                    
            except Exception as e:
                logger.error(f"Client error: {e}")
            finally:
                writer.close()
                await writer.wait_closed()
                logger.info(f"Client disconnected from {addr}")
        
        server = await asyncio.start_server(
            handle_client, host, port
        )
        
        async with server:
            await server.serve_forever()


def mcp_tool(params: Dict[str, Any] = None):
    """Decorator to mark methods as MCP tools"""
    def decorator(func):
        func._mcp_params = params or {}
        return func
    return decorator
