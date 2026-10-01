import sys, os, io, csv, zipfile, json, time, pickle
from pathlib import Path
sys.path.insert(0, os.path.expandvars("$EVITEA_REPO/src")); sys.path.insert(0, os.path.expandvars("$EVITEA_REPO"))
os.environ["PYTHONDONTWRITEBYTECODE"]="1"; sys.dont_write_bytecode=True
AKA=Path(os.path.expandvars("$EVITEA_DATA/AKA017/route_b_r9e_validation_full55x/shadow_v1/reports"))
ASSET=Path(os.path.expandvars("$EVITEA_DATA/TEMP2_ISO1"))
TE_FASTA=ASSET/"TE consensus library"/"transposon_sequence_set.fa"
FAMMAP=ASSET/"fbte_to_flybase_family.tsv"
WORK=Path(os.path.expanduser("~/scratch/g1"))
ZIPS={"g1a":"G1A_native_results/TE_LOCATOR_ROUTE_B_G1A_RESULTS_FOR_CHATGPT.zip","g1b":"G1B_sparse_panel_routing/TE_LOCATOR_ROUTE_B_G1B_RESULTS_FOR_CHATGPT.zip","g1c":"G1C_coverage_rescue/TE_LOCATOR_ROUTE_B_G1C_RESULTS_FOR_CHATGPT.zip","g1d":"G1D_optimized_rescue_routing/TE_LOCATOR_ROUTE_B_G1D_RESULTS_FOR_CHATGPT.zip","g1e":"G1E_pair_aware_rescue/TE_LOCATOR_ROUTE_B_G1E_RESULTS_FOR_CHATGPT.zip","g1f":"G1F_host_mate_locus_gate/TE_LOCATOR_ROUTE_B_G1F_RESULTS_FOR_CHATGPT.zip","g1g":"G1G_reference_resident_compatibility/TE_LOCATOR_ROUTE_B_G1G_RESULTS_FOR_CHATGPT.zip","g1hr":"G1HR_candidate_recall_repair/TE_LOCATOR_ROUTE_B_G1HR_RESULTS_FOR_CHATGPT.zip"}
def ztsv(stage, member):
    with zipfile.ZipFile(AKA/ZIPS[stage]) as z, z.open(member) as raw:
        return list(csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8", newline=""), delimiter="\t"))
def znpy(stage, member):
    import numpy as np
    with zipfile.ZipFile(AKA/ZIPS[stage]) as z: return np.load(io.BytesIO(z.read(member)), allow_pickle=False)
def zjson(stage, member):
    with zipfile.ZipFile(AKA/ZIPS[stage]) as z: return json.loads(z.read(member))
def cached(name, fn):
    p=WORK/(name+".pkl")
    if p.exists(): return pickle.loads(p.read_bytes())
    v=fn(); p.write_bytes(pickle.dumps(v)); return v
