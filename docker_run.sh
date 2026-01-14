#!/bin/bash

# Mount .env file if it exists
ENV_MOUNT=""
if [ -f ".env" ]; then
    ENV_MOUNT="-v $(pwd)/.env:/app/.env"
fi

# Run with data volume mounted for caching
docker run -it \
    -p 8501:8501 \
    -v "$(pwd)/data:/app/data" \
    $ENV_MOUNT \
    natalie23gill/geneset-translator:1.0.0 \
    streamlit run app.py --server.address 0.0.0.0
