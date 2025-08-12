"""
This is an example dataset.py file that can be loaded dynamically by bonito
"""

from functools import partial
from pathlib import Path

import numpy as np
from torch.utils.data import RandomSampler

from bonito.data import ChunkDataSet

class ChunkLoader:

    def __init__(self, train_data, valid_data, **kwargs):
        self.train_data = train_data
        self.valid_data = valid_data

    def train_loader_kwargs(self, **kwargs):
        train_ds = ChunkDataSet(*self.train_data)
        return {
            "dataset": train_ds,
            "sampler": RandomSampler(train_ds, num_samples=kwargs["chunks"]),
        }

    def valid_loader_kwargs(self, **kwargs):
        valid_ds = ChunkDataSet(*self.valid_data)
        return {
            "dataset": valid_ds,
            "shuffle": False,
        }


def load_chunks(input_folder):
    chunks = np.load(input_folder / "chunks.npy").astype(np.float32)
    refs = np.load(input_folder / "references.npy").astype(np.int64)
    ref_lens = np.load(input_folder / "reference_lengths.npy").astype(np.int64)
    return chunks, refs, ref_lens


chunks, refs, lens = load_chunks(Path("/data/ctc_output"))

# As an example, we take the first 1000 chunks for training and the last 100 for validation
# In practice more data will be required! 
train_data = chunks[:1000], refs[:1000], lens[:1000]
valid_data = chunks[100:], refs[100:], lens[100:]

Loader = partial(ChunkLoader, train_data, valid_data)
