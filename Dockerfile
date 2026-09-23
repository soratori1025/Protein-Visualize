FROM python:3.11-slim

# Cài đặt DSSP từ repo chính thức của Debian (chuẩn nhất, không bị lỗi C++ ABI) và curl/bzip2
RUN apt-get update && apt-get install -y \
    dssp curl bzip2 \
    && rm -rf /var/lib/apt/lists/*

# Tải file từ điển nén và giải nén thành components.cif tĩnh tại /var/cache/libcifpp (đúng như DSSP yêu cầu)
RUN mkdir -p /var/cache/libcifpp && \
    curl -s -o /var/cache/libcifpp/components.cif.gz https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz && \
    gunzip /var/cache/libcifpp/components.cif.gz && \
    curl -s -o /var/cache/libcifpp/mmcif_pdbx.dic https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic

# Cài đặt STRIDE qua Micromamba (vì apt không có stride) và đưa vào PATH
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /usr/local bin/micromamba && \
    /usr/local/bin/micromamba create -y -p /opt/conda -c conda-forge -c bioconda stride && \
    ln -s /opt/conda/bin/stride /usr/bin/stride

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
