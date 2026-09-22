FROM python:3.11-slim

# Cài đặt DSSP
RUN apt-get update && apt-get install -y \
    dssp \
    && rm -rf /var/lib/apt/lists/*

# Thiết lập thư mục gốc của repo
WORKDIR /workspace

# Chỉ copy requirements trước để tận dụng cache của Docker
COPY backend/requirements.txt ./backend/
RUN pip install --no-cache-dir -r backend/requirements.txt

# Copy mã nguồn backend và protein_engine (đây là mấu chốt để fix lỗi ModuleNotFoundError)
COPY backend ./backend
COPY protein_engine ./protein_engine

# Chuyển thư mục làm việc vào backend để giống môi trường local của bạn
WORKDIR /workspace/backend

# Thêm workspace vào PYTHONPATH để Python nhận diện được gói protein_engine
ENV PYTHONPATH=/workspace

# Khởi chạy Uvicorn
CMD uvicorn app.main:app --host 0.0.0.0 --port ${PORT:-8000}
