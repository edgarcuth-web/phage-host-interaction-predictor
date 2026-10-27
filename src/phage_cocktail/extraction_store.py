"""Small immutable object store for extraction checkpoints, not raw tensors."""
from pathlib import Path, PurePosixPath
import os
import tempfile
from urllib.parse import urlsplit


class StoreError(RuntimeError):
    pass


def object_key(key):
    if (not isinstance(key, str) or not key or key.startswith("/") or "\\" in key
            or any(p in ("", ".", "..") for p in key.split("/"))):
        raise ValueError("invalid checkpoint object key")
    return str(PurePosixPath(key))


class LocalStore:
    """Publish complete files by atomic hard-link creation on a local filesystem."""
    def __init__(self, root):
        self.root = Path(root).resolve()

    def _path(self, key):
        path = (self.root / object_key(key)).resolve()
        if not path.is_relative_to(self.root):
            raise ValueError("checkpoint path escapes its root")
        return path

    def get(self, key):
        try:
            return self._path(key).read_bytes()
        except FileNotFoundError:
            return None

    def create(self, key, payload):
        path = self._path(key)
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
                temporary = Path(stream.name)
                stream.write(payload)
                stream.flush()
                os.fsync(stream.fileno())
            try:
                os.link(temporary, path)
            except FileExistsError:
                return False
            directory = os.open(path.parent, os.O_RDONLY)
            try:
                os.fsync(directory)
            finally:
                os.close(directory)
            return True
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)


class S3Store:
    """Require real conditional PUT semantics; never fall back to overwriting."""
    def __init__(self, bucket, prefix, *, client=None, endpoint_url=None):
        if not bucket or not prefix:
            raise ValueError("use a dedicated non-empty S3 run prefix")
        self.bucket, self.prefix = bucket, object_key(prefix.strip("/"))
        if endpoint_url is not None:
            endpoint = urlsplit(endpoint_url)
            if (endpoint.scheme != "https" or not endpoint.netloc or endpoint.username
                    or endpoint.password or endpoint.query or endpoint.fragment):
                raise ValueError("S3 endpoint must be an HTTPS URL without credentials or query")
        if client is None:
            try:
                import boto3
                from botocore.config import Config
            except ImportError:
                raise ValueError("S3 requires pip install -e '.[cloud]'") from None
            client = boto3.client("s3", endpoint_url=endpoint_url, config=Config(
                signature_version="s3v4", connect_timeout=10, read_timeout=60,
                retries={"total_max_attempts": 1}))
        self.client = client

    def check_conditional_writes(self):
        key, payload = "conditional_write_probe", b"immutable-checkpoint-probe-v1"
        self.create(key, payload)
        if self.create(key, b"overwrite-must-fail") or self.get(key) != payload:
            raise StoreError("object store does not provide required conditional-create semantics")

    def _key(self, key):
        return self.prefix + "/" + object_key(key)

    def get(self, key):
        try:
            response = self.client.get_object(Bucket=self.bucket, Key=self._key(key))
            with response["Body"] as stream:
                payload = stream.read(64 * 1024 * 1024 + 1)
            if len(payload) > 64 * 1024 * 1024:
                raise StoreError("checkpoint object exceeds 64 MiB")
            return payload
        except Exception as exc:
            code = getattr(exc, "response", {}).get("Error", {}).get("Code")
            if code in ("NoSuchKey", "404"):
                return None
            raise StoreError("S3 checkpoint read failed") from None

    def create(self, key, payload):
        try:
            self.client.put_object(Bucket=self.bucket, Key=self._key(key), Body=payload,
                                   IfNoneMatch="*", ContentType="application/octet-stream")
            return True
        except Exception as exc:
            code = getattr(exc, "response", {}).get("Error", {}).get("Code")
            if code in ("PreconditionFailed", "412"):
                return False
            # Timeouts/409s are ambiguous. Stop, and inspect the persisted object.
            raise StoreError("S3 checkpoint create failed; no application retry") from None


def open_store(location, *, endpoint_url=None):
    if location.startswith("s3://"):
        parsed = urlsplit(location)
        if parsed.query or parsed.fragment or parsed.username or parsed.password:
            raise ValueError("S3 location must contain only bucket and run prefix")
        return S3Store(parsed.netloc, parsed.path.strip("/"), endpoint_url=endpoint_url)
    if "://" in location or endpoint_url is not None:
        raise ValueError("use a local directory or s3://bucket/run-prefix")
    return LocalStore(location)


def put_identical(store, key, payload):
    """Idempotent derived/imported artifacts; live request claims use create directly."""
    if not store.create(key, payload) and store.get(key) != payload:
        raise ValueError(f"checkpoint object differs: {key}")
