FROM python:3.11-slim
RUN apt-get update && apt-get install -y --no-install-recommends libgl1 && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
COPY data/reference ./data/reference
COPY data/manifests ./data/manifests
RUN pip install --no-cache-dir -e ".[slides]" --extra-index-url https://download.pytorch.org/whl/cpu
ENTRYPOINT ["hrdscope"]
CMD ["--help"]
