import gc
import os
import numpy as np

import MDAnalysis as mda

from tqdm import tqdm


class BatchLoader:
    def __init__(self,
                 reference,
                 trj_list,
                 trajectory_stride,
                 batch_size,
                 dt_ns,
                 transforms=None,
                 topologyattrs=None,
                 pattern="all",
                 ):
        self.reference = reference
        self.trj_list = trj_list
        self.trajectory_stride = trajectory_stride
        self.batch_size = batch_size
        self.dt_ns = dt_ns

        self.transforms = transforms
        self.topologyattrs = topologyattrs
        self.pattern = pattern

    def load(self):
        for trj_split_index in tqdm(range(0, len(self.trj_list), self.batch_size)):
            batch_index = trj_split_index // self.batch_size

            batch_u = mda.Universe(self.reference._topology,
                                   self.trj_list[trj_split_index:trj_split_index + self.batch_size],
                                   in_memory=True,
                                   in_memory_step=self.trajectory_stride,
                                   topology_format="PDB",
                                   )
            if self.topologyattrs:
                for key, value in self.topologyattrs.items():
                    batch_u.add_TopologyAttr(key, value)

            if self.transforms:
                transforms = [transform(batch_u.atoms) for transform in self.transforms]
                batch_u.trajectory.add_transformations(*transforms)

            batch_ag = batch_u.select_atoms(self.pattern)
            indexes = np.linspace(0, batch_ag.universe.trajectory.n_frames * self.trajectory_stride,
                                  batch_ag.universe.trajectory.n_frames,
                                  endpoint=False)
            indexes += (indexes[-1] + self.trajectory_stride) * (batch_index)
            batch_time_ns = indexes * self.dt_ns

            yield batch_ag, batch_time_ns


class BatcHist3dWriter:

    def __init__(self, output_directory=".", header=None):
        self.output_directory = output_directory
        self.header = header
        self.results = None

    def save(self, batch_time_ns, batch_result):
        os.makedirs(self.output_directory, exist_ok=True)
        if self.results is None:
            self.results = batch_result.results.density.grid
            np.save(os.path.join(self.output_directory, "edges.npy"), batch_result.results.density.edges)
        else:
            self.results += batch_result.results.density.grid
        np.save(os.path.join(self.output_directory, "hist.npy"), self.results)


class BatchAnalyzer:
    def __init__(self,
                 batchloader,
                 analyzer,
                 writer):
        self.batchloader = batchloader
        self.analyzer = analyzer
        self.writer = writer

    def analyse(self):
        for batch_u, batch_time_ns in self.batchloader.load():
            result = self.analyzer(batch_u).run()
            self.writer.save(batch_time_ns, result)
            gc.collect()
