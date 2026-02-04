import numpy as np

from MDAnalysis.analysis.density import DensityAnalysis


class AnalyzerWrapper:
    def __init__(self, analysis, *args, **kwargs):
        self.analysis = analysis

        self.__dict__.update(kwargs)
        self.attr_names = list(kwargs.keys())

        for arg in args:
            setattr(self, arg, arg)
            self.attr_names.append(arg)

    def __call__(self, ag):
        attrs = [getattr(self, arg) for arg in self.attr_names]
        return self.analysis(ag, *attrs)


class Hist3dAnalysis(DensityAnalysis):

    def _conclude(self):
        super(Hist3dAnalysis, self)._conclude()

        dedges = [np.diff(edge) for edge in self.results.density.edges]
        D = len(self.results.density.edges)
        for i in range(D):
            shape = np.ones(D, int)
            shape[i] = len(dedges[i])
            self.results.density.grid *= dedges[i].reshape(shape)

        self.results.density.grid *= self.n_frames
