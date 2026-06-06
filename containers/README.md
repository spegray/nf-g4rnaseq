# Container image

The R/analysis processes run in a single image built from `Dockerfile`
(+ `environment.yml`). The standard tools (STAR, fastp, samtools, featureCounts,
umi_tools, MultiQC) use their own biocontainers per process, so only this image
needs publishing. The environment is verified to solve on **linux-64**
(361 packages) — the image builds cleanly.

## Automatic build & push (recommended — no local Docker needed)

A version tag triggers `.github/workflows/build-container.yml`, which builds for
linux/amd64 and pushes to the GitHub Container Registry:

```bash
git tag v0.1.0 && git push origin v0.1.0
# -> ghcr.io/<your-github-owner>/nf-g4rnaseq-r:0.1.0  (and :latest)
```

Then make the package **public** (GitHub repo → Packages → package → settings)
and point the pipeline at it — either edit `nextflow.config`:

```groovy
params.r_container = "ghcr.io/<your-github-owner>/nf-g4rnaseq-r:0.1.0"
```

or pass it per run: `--r_container ghcr.io/<owner>/nf-g4rnaseq-r:0.1.0`.

## Manual build (if you have Docker locally)

```bash
docker buildx build --platform linux/amd64 \
  -t ghcr.io/<owner>/nf-g4rnaseq-r:0.1.0 -f Dockerfile .
docker push ghcr.io/<owner>/nf-g4rnaseq-r:0.1.0
```

## Notes
- `linux/amd64` is the nf-core / HPC / cloud default. An arm64 build is optional
  (a few Bioconductor builds differ); add `linux/arm64` to `platforms:` if needed.
- `-profile conda` uses this same `environment.yml` directly (no image required),
  which is also what the CI end-to-end test uses.
