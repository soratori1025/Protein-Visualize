from pathlib import Path


class ProteinService:
    @staticmethod
    def ensure_upload_dir() -> Path:
        upload_dir = Path("data/uploads")
        upload_dir.mkdir(parents=True, exist_ok=True)
        return upload_dir

    @staticmethod
    def save_uploaded_file(file_content: bytes, filename: str) -> Path:
        upload_dir = ProteinService.ensure_upload_dir()
        path = upload_dir / filename
        path.write_bytes(file_content)
        return path