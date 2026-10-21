"""
preprocess.py -- Cirro dataset preprocessing hook for the CRAM-to-BAM
conversion pipeline. Mirrors the MuTect1 pipeline's preprocess.py pattern:
injects params.cram_runs directly, no samplesheet CSV involved.

Every .cram / .cram.crai in the dataset is considered. If the user supplies
the optional 'cram_folder' parameter, only files under that folder (relative
to the dataset's data/ root) are used; otherwise the whole dataset is searched.
"""
import json
from cirro.helpers.preprocess_dataset import PreprocessDataset


def relative_path(uri):
    """Path of a file relative to the dataset's data/ root."""
    marker = "/data/"
    idx = uri.find(marker)
    return uri[idx + len(marker):] if idx != -1 else uri


def extract_crams(ds, cram_folder=""):
    df = ds.files.copy()
    df["file"] = df["file"].astype(str)

    df = df[df["file"].str.endswith(".cram") | df["file"].str.endswith(".cram.crai")]

    folder = (cram_folder or "").strip().strip("/")
    if folder:
        prefix = folder + "/"
        df = df[df["file"].map(lambda f: relative_path(f).startswith(prefix))]

    cram_map = {}

    for sample, group in df.groupby("sample"):
        cram = ""
        crai = ""

        for f in group["file"]:
            if f.endswith(".cram") and not f.endswith(".cram.crai"):
                if cram:
                    raise ValueError(
                        f"Multiple CRAMs found for sample {sample!r}: "
                        f"{cram!r} and {f!r} -- expected exactly one. "
                        "Use the 'cram_folder' parameter to narrow the search."
                    )
                cram = f
            elif f.endswith(".cram.crai"):
                crai = f

        if cram:
            cram_map[str(sample)] = {"cram": cram, "crai": crai}

    if not cram_map:
        where = f"under folder {folder!r}" if folder else "anywhere"
        raise ValueError(f"No CRAMs found {where} in this dataset")

    return cram_map


def main():
    ds = PreprocessDataset.from_running()

    print("=== ds.files preview ===")
    print(ds.files.head(20).to_string(index=False))

    cram_folder = ds.params.get("cram_folder", "")
    print(f"cram_folder: {cram_folder!r}")
    cram_map = extract_crams(ds, cram_folder)

    cram_runs = []
    for sample, files in cram_map.items():
        if not files["crai"]:
            raise ValueError(f"Sample {sample!r} has a CRAM but no matching .crai index")

        cram_runs.append({
            "sample_id": sample,
            "cram": files["cram"],
            "crai": files["crai"],
        })

    ds.add_param("cram_runs", cram_runs)

    print("\nFinal parameters:")
    print(json.dumps(ds.params, indent=2, default=str))


if __name__ == "__main__":
    main()