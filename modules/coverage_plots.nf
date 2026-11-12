process COVERAGE_PLOTS {
    label 'process_low'
    conda "${projectDir}/envs/python_analysis_env.yaml"
    errorStrategy "ignore"
    publishDir "${params.outdir}/08_coverage/plots", mode: params.publish_dir_mode

    input:
    path depth_files

    output:
    path "*.png", emit: plots
    path "versions.yml", emit: versions

    script:
    """
    python3 - << 'PY'
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

files = list(Path(".").glob("*.depth"))
if not files:
    Path("no_data.png").write_bytes(b"")
else:
    # per-sample plots
    for f in files:
        sample = f.stem
        pos, depth = [], []
        with f.open() as fh:
            for line in fh:
                p = line.strip().split()
                if len(p) >= 3:
                    pos.append(int(p[1])); depth.append(int(p[2]))
        if not pos:
            continue
        fig, ax = plt.subplots(figsize=(12,3))
        ax.fill_between(pos, depth, alpha=0.6)
        ax.set_xlabel("Genome position"); ax.set_ylabel("Depth")
        ax.set_title(f"Coverage: {sample}")
        ax.axhline(10, color="red", ls="--", lw=0.8, label="10x")
        ax.legend(loc="upper right")
        fig.tight_layout()
        fig.savefig(f"{sample}_coverage.png", dpi=120)
        plt.close()

    # multi-sample heatmap (downsample)
    if len(files) > 1:
        data, labels = [], []
        min_len = None
        for f in files:
            d = []
            with f.open() as fh:
                for line in fh:
                    p = line.strip().split()
                    if len(p) >= 3: d.append(int(p[2]))
            if not d: continue
            if min_len is None or len(d) < min_len: min_len = len(d)
            data.append(d); labels.append(f.stem)
        if data and min_len:
            step = max(1, min_len // 500)
            mat = np.array([d[:min_len:step] for d in data])
            fig, ax = plt.subplots(figsize=(14, max(3, len(labels)*0.4)))
            im = ax.imshow(np.log1p(mat), aspect="auto", interpolation="nearest")
            ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels, fontsize=7)
            ax.set_xlabel("Genome position (downsampled)")
            ax.set_title("Coverage heatmap (log1p depth)")
            fig.colorbar(im, ax=ax, fraction=0.02)
            fig.tight_layout()
            fig.savefig("coverage_heatmap.png", dpi=120)
            plt.close()
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
