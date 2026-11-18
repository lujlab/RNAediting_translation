from pathlib import Path
import csv, subprocess, concurrent.futures

p = Path(__file__).resolve().parent
jobs = sorted(
    {
        (r["species"], r["dataset"], r["library"])
        for r in csv.DictReader(open(p / "inputs.tsv"), delimiter="\t")
        if r["species"] in ("human", "mouse")
    }
)


def run(job):
    prefix = ".".join(job)
    with open(p / (prefix + ".all_CDS.log"), "w") as log:
        subprocess.run(
            [str(p / "count_all_CDS"), str(p), *job], stdout=log, stderr=log, check=True
        )
    print("Finished", prefix, flush=True)


with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    list(pool.map(run, jobs))
