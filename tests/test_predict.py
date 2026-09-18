import json
import h5py
import numpy as np
import torch

from hrdscope.mil import GatedAttentionMIL
from hrdscope.predict import predict_features


def test_predict_features_ensemble(tmp_path):
    run = tmp_path / "run"; run.mkdir()
    for seed in (0, 1):
        torch.save(GatedAttentionMIL(8).state_dict(), run / f"abmil_seed{seed}_fold0.pt")
    with h5py.File(tmp_path / "s.h5", "w") as h5:
        h5.create_dataset("features", data=np.random.default_rng(0).normal(size=(30, 8)).astype(np.float16))
        h5.create_dataset("coords", data=np.zeros((30, 2), dtype=np.int32))
        h5.attrs.update({"slide": "s.svs", "backbone": "test", "tile_size0": 448})
    res = predict_features(tmp_path / "s.h5", run)
    assert res["n_models"] == 2 and 0 <= res["probability"]["ge42"]["mean"] <= 1
    assert len(res["attention"]["weights"]) == 30
    json.dumps(res)
