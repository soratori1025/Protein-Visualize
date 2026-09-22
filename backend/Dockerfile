FROM python:3.11-slim

# Cài đặt DSSP để thuật toán dssp_ss và dssp_slab hoạt động
RUN apt-get update && apt-get install -y \
    dssp \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy toàn bộ code backend vào container
COPY . .

# Render sẽ tự động gán biến môi trường PORT
CMD uvicorn app.main:app --host 0.0.0.0 --port ${PORT:-8000}
