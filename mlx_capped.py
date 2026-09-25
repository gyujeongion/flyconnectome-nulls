"""Run an engine with MLX's buffer cache capped (MLX_CACHE_GB, default 2). The cache only keeps freed buffers for
reuse; capping it changes memory use, not any computation. Usage: python mlx_capped.py evo_bsham.py <engine args>"""
import os, runpy, sys
import mlx.core as mx
mx.set_cache_limit(int(float(os.environ.get("MLX_CACHE_GB", "2")) * 2**30))
engine = sys.argv[1]; sys.argv = sys.argv[1:]
try:
    runpy.run_path(engine, run_name="__main__")
finally:
    print(f"mlx memory: peak {mx.get_peak_memory() / 2**30:.1f} GB, cache limit {os.environ.get('MLX_CACHE_GB', '2')} GB", flush=True)
