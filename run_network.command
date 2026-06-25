#!/bin/zsh
cd "$(dirname "$0")"
HOST=0.0.0.0 PORT=8788 /Users/zhangyu/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 server.py
