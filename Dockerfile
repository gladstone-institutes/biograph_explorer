FROM python:3.11-slim

WORKDIR /app

# Install Poetry
RUN pip install poetry
# Copy application code
COPY . .

# Install dependencies (no dev dependencies, no virtualenv in container)
RUN poetry config virtualenvs.create false && \
    poetry install --no-interaction --no-ansi --only main


CMD ["streamlit", "run", "app.py"]
