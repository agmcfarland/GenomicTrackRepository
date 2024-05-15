
import warnings
from os.path import join as pjoin
from os.path import exists

class ProjectPaths:

	def __init__(self, base_directory : str):

		self.base_directory = base_directory

		self._generate_path_skeleton()

	def _generate_path_skeleton(self) -> dict:

		self.paths = {
			'base' : self.base_directory,

			'data' :  pjoin(self.base_directory, 'data'),
			'data_raw' : pjoin(self.base_directory, 'data', 'raw'),
			'data_processed' : pjoin(self.base_directory, 'data', 'processed'),
			'data_interim' : pjoin(self.base_directory, 'data', 'interim'),
			'data_external' : pjoin(self.base_directory, 'data', 'external'),

			'docs' : pjoin(self.base_directory, 'docs'),

			'environment' : pjoin(self.base_directory, 'environment'),

			'exploratory' : pjoin(self.base_directory, 'exploratory'),

			'models' : pjoin(self.base_directory, 'models'),

			'notebooks' : pjoin(self.base_directory, 'notebooks'),

			'references' : pjoin(self.base_directory, 'references'),

			'reports' : pjoin(self.base_directory, 'reports'),

			'src' : pjoin(self.base_directory, 'src'),

			'tests' : pjoin(self.base_directory, 'tests')

		}

	def add_path(self, key : str, path : str) -> str:

		self.paths[key] = path

		if not exists(path):
			warnings.warn(f'does not exist: {path}', UserWarning)


