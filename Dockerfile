FROM python:3.11-slim

WORKDIR /app

# Install uv
RUN pip install uv
# Copy application code
COPY . .

# Install the project and its runtime dependencies into the system environment
# (no dev dependencies, no virtualenv in container)
RUN uv pip install --system --no-cache .


CMD ["streamlit", "run", "app.py"]
