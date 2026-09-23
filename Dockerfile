FROM python:3.11-slim

# Cài đặt curl và bzip2 để tải Micromamba
RUN apt-get update && apt-get install -y \
    curl bzip2 \
    && rm -rf /var/lib/apt/lists/*

# Thiết lập thư mục gốc của repo
WORKDIR /workspace

# Cài đặt Micromamba, DSSP 3.0.0 (từ salilab) và STRIDE (từ bioconda)
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /usr/local bin/micromamba && \
    /usr/local/bin/micromamba create -y -p /workspace/tools/conda-env -c conda-forge -c salilab -c bioconda dssp=3.0.0 stride && \
    mkdir -p /workspace/tools/bin && \
    echo '#!/usr/bin/env bash\nexport LD_LIBRARY_PATH="/workspace/tools/conda-env/lib:$LD_LIBRARY_PATH"\nexec "/workspace/tools/conda-env/bin/mkdssp" "$@"' > /workspace/tools/bin/dssp && \
    echo '#!/usr/bin/env bash\nexport LD_LIBRARY_PATH="/workspace/tools/conda-env/lib:$LD_LIBRARY_PATH"\nexec "/workspace/tools/conda-env/bin/mkdssp" "$@"' > /workspace/tools/bin/mkdssp && \
    echo '#!/usr/bin/env bash\nexport LD_LIBRARY_PATH="/workspace/tools/conda-env/lib:$LD_LIBRARY_PATH"\nexec "/workspace/tools/conda-env/bin/stride" "$@"' > /workspace/tools/bin/stride && \
    chmod +x /workspace/tools/bin/*

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
