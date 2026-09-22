import os
import sys
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

TOOLS_BIN = ROOT / "tools" / "bin"
TOOLS_SHARE = ROOT / "tools" / "share" / "libcifpp"
if TOOLS_BIN.is_dir() and str(TOOLS_BIN) not in os.environ.get("PATH", ""):
    os.environ["PATH"] = f"{TOOLS_BIN}{os.pathsep}{os.environ.get('PATH', '')}"
if TOOLS_SHARE.is_dir() and "LIBCIFPP_DATA_DIR" not in os.environ:
    os.environ["LIBCIFPP_DATA_DIR"] = str(TOOLS_SHARE)

from app.api.structure import router as structure_router
from app.api.secondary_structure import router as secondary_structure_router
from app.api.analysis import router as analysis_router

app = FastAPI(
    title="ProteinLab API",
    version="0.1.0",
    description="MVP backend for protein structure upload and analysis.",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], # Allow any frontend (like Vercel) to call this API
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/")
def read_root():
    return {"message": "ProteinLab API is running on Render!"}

@app.get("/api/health")
def health() -> dict:
    return {
        "status": "ok",
        "service": "ProteinLab API",
        "version": "0.1.0",
    }


app.include_router(structure_router)
app.include_router(secondary_structure_router)
app.include_router(analysis_router)
