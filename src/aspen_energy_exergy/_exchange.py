"""Case identity and atomic publication for Aspen exchange documents."""

import hashlib
import json
import os
from pathlib import Path
import re
from tempfile import NamedTemporaryFile


def digest(value):
    data = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def case_hashes(model):
    meta = model["metadata"]
    values = [meta.get(k) for k in ("backup_sha256", "source_sha256", "archive_sha256")]
    hashes = [h.lower() for h in values if isinstance(h, str) and re.fullmatch("[a-fA-F0-9]{64}", h)]
    if hashes:
        return hashes
    inventory = {key: {row["id"]: row for row in model[key]} for key in ("components", "streams")}
    return [digest(inventory)]


def publish_template(destination, data):
    """Create the complete file once or reuse an identical copy, preserving filled returns."""
    path = Path(destination)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with NamedTemporaryFile(dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(data)
        if os.name == "nt":
            temporary.rename(path)
        else:
            os.link(temporary, path)
    except FileExistsError:
        if path.read_bytes() != data:
            raise FileExistsError(f"A different file already exists at {path}. Choose a new template filename.") from None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return path
