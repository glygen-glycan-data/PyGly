FROM node:20-slim

RUN apt-get update && apt-get install -y \
    curl \
    git \
    procps \
    openssh-client \
    python3 \
    && rm -rf /var/lib/apt/lists/*

RUN npm install -g @anthropic-ai/claude-code@2.1.112
USER 1000
WORKDIR /workspace
