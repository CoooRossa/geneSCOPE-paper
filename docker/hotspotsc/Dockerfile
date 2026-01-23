FROM python:3.11-slim

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    HOME=/tmp \
    NUMBA_CACHE_DIR=/tmp

# Minimal runtime libs (OpenMP often needed by numpy/scipy wheels)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgomp1 \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/app
COPY requirements.txt /opt/app/requirements.txt
RUN pip install --no-cache-dir -r /opt/app/requirements.txt

# Copy runnable script and entrypoint
COPY app/run_hotspot_xenium.py /opt/app/run_hotspot_xenium.py
COPY bin/entrypoint.sh /opt/app/entrypoint.sh
RUN chmod +x /opt/app/entrypoint.sh

WORKDIR /work
ENTRYPOINT ["/opt/app/entrypoint.sh"]
