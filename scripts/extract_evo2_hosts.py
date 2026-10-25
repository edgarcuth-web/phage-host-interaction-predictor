"""Freeze a local 40B deployment, then extract pilot/full host caches offline."""
import argparse
from hashlib import sha256
from importlib.metadata import version
import inspect
from pathlib import Path
import pkgutil

from phage_cocktail.artifacts import atomic_write_json, read_json_object, sha256_file
from phage_cocktail.evo2_host import PROTOCOL, extract_cache, validate_deployment


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True, help="Existing local merged checkpoint; never downloaded")
    parser.add_argument("--deployment-lock", type=Path, required=True)
    parser.add_argument("--freeze-deployment", action="store_true", help="Write a new lock and exit before extracting")
    parser.add_argument("--pilot", action="store_true", help="First three sorted hosts, without assay labels")
    parser.add_argument("--data-root", type=Path, default=Path("data/raw/genophi_host_study"))
    parser.add_argument("--outdir", type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    # GPU dependencies are optional and stay off the benchmark/test dependency path.
    import torch
    import transformer_engine
    from evo2 import Evo2
    from vortex.model.model import StripedHyena

    if not torch.cuda.is_available():
        raise ValueError("40B extraction requires the deployed CUDA environment")
    spec = read_json_object(repo / PROTOCOL)["embedding"]
    config = pkgutil.get_data("evo2", "configs/evo2-40b-1m.yml")
    identity = {
        "model_name": spec["model_name"], "checkpoint_sha256": sha256_file(args.checkpoint),
        "config_sha256": sha256(config).hexdigest(),
        "evo2_module_sha256": sha256_file(inspect.getfile(Evo2)),
        "vortex_module_sha256": sha256_file(inspect.getfile(StripedHyena)),
        "evo2_version": version("evo2"), "vortex_version": version("vtx"),
        "torch_version": torch.__version__, "cuda_version": torch.version.cuda,
        "transformer_engine_version": transformer_engine.__version__,
        "hardware": "|".join(torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())),
        "precision": "Upstream checkpoint/config precision with FP8 input projections; no override or autocast; FP32 token means",
    }
    validate_deployment(identity)
    if args.freeze_deployment:
        if args.deployment_lock.exists():
            raise ValueError("deployment lock must be new")
        atomic_write_json(identity, args.deployment_lock)
        return
    if read_json_object(args.deployment_lock) != identity:
        raise ValueError("runtime differs from frozen deployment; amend before extracting")
    if args.outdir is None:
        raise ValueError("extraction requires --outdir")
    model = Evo2(spec["model_name"], local_path=str(args.checkpoint), use_kernels=False)
    model.model.eval()
    model.model.get_submodule(spec["layer"])  # Fail early; never choose another layer.
    for i in range(torch.cuda.device_count()):
        torch.cuda.reset_peak_memory_stats(i)

    def encode(sequence):
        tokens = model.tokenizer.tokenize(sequence)
        if len(tokens) != len(sequence):
            raise ValueError("protocol requires one token per input base, without special tokens")
        ids = torch.tensor(tokens, dtype=torch.int, device="cuda:0").unsqueeze(0)
        with torch.inference_mode():
            _, embeddings = model(ids, return_embeddings=True, layer_names=[spec["layer"]])
        values = embeddings[spec["layer"]]
        if tuple(values.shape) != (1, len(sequence), spec["dimension"]):
            raise ValueError("40B layer output shape differs from protocol")
        return values[0].float().cpu().numpy()

    def runtime_metrics():
        return {"cuda_peak_memory_bytes": [
            {"device": i, "allocated": torch.cuda.max_memory_allocated(i),
             "reserved": torch.cuda.max_memory_reserved(i)} for i in range(torch.cuda.device_count())]}

    extract_cache(repo, args.data_root, args.outdir, encode, identity,
                  pilot=args.pilot, runtime_metrics=runtime_metrics)


if __name__ == "__main__":
    main()
