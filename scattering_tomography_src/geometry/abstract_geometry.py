from abc import ABC, abstractmethod

class AbstractGeometry(ABC):
    """
    Inteface for the geometry class, which contains all the elements of the geometry - detectors, sources and volume.
    """
    @abstractmethod
    def __init__(self):
        self._detectors = None
        self._sources = None
        self._volume = None


