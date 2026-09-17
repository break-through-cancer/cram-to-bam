"""
preprocess.py -- Cirro dataset preprocessing hook for the CRAM-to-BAM
conversion pipeline. Mirrors the MuTect1 pipeline's preprocess.py pattern:
injects params.cram_runs directly, no samplesheet CSV involved.

Only CRAMs living under a 'recalibrated/' folder are considered -- datasets
may contain CRAMs from earlier, non-recalibrated pipeline stages too, and
those must be excluded.
"""
import json
from cirro.helpers.preprocess_dataset import PreprocessDataset

RECAL_FOLDER_PATTERN = r"/recalibrated/"


def extract_crams(ds):
    df = ds.files.copy()
    df["file"] = df["file"].astype(str)

    # Restrict to the recalibrated folder BEFORE anything else -- this is
    # the filter that keeps this dataset's output limited to the CRAMs we
    # actually want, even if other pipeline stages' CRAMs are also present.
    df = df[df["file"].str.contains(RECAL_FOLDER_PATTERN, regex=True)]

    df = df[df["file"].str.endswith(".cram") | df["file"].str.endswith(".cram.crai")]

    cram_map = {}

    for sample, group in df.groupby("sample"):
        cram = ""
        crai = ""

        for f in group["file"]:
            if f.endswith(".cram") and not f.endswith(".cram.crai"):
                if cram:
                    raise ValueError(
                        f"Multiple recalibrated CRAMs found for sample {sample!r}: "
                        f"{cram!r} and {f!r} -- expected exactly one."
                    )
                cram = f
            elif f.endswith(".cram.crai"):
                crai = f

        if cram:
            cram_map[str(sample)] = {"cram": cram, "crai": crai}

    if not cram_map:
        raise ValueError(
            f"No CRAMs found under a folder matching {RECAL_FOLDER_PATTERN!r} in this dataset"
        )

    return cram_map


def main():
    ds = PreprocessDataset.from_running()

    print("=== ds.files preview ===")
    print(ds.files.head(20).to_string(index=False))

    cram_map = extract_crams(ds)

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