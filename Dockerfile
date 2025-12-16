FROM ubuntu:25.10

# Install dependencies
RUN apt-get update && apt-get install -y \
    python3-dev \
    gcc \
    python3-venv

RUN apt-get clean
RUN rm -rf /var/lib/apt/lists/*

COPY requirements.txt .

# Install packages with CodeArtifact as additional index
RUN python3 -m venv artemis
RUN artemis/bin/pip3 install --no-cache-dir \
    -r requirements.txt

RUN find . -name "__pycache__" -type d -exec rm -rf {} +
RUN rm -rf *.egg-info

# Activate the virtual environment by default
ENV PATH="/artemis/bin:$PATH"