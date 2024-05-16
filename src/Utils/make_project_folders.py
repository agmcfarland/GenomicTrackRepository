
import os
from os.path import join as pjoin
from src.ProjectPaths import ProjectPaths

# Base path
project_paths = ProjectPaths.ProjectPaths('/data/GenomicTrackRepository')


os.makedirs(pjoin(project_paths.paths['notebooks'], 'hg38'), exist_ok = True)

# iGUIDE software location
project_paths.add_path(key = 'iGUIDE_path', path = '/data/iGUIDE')

# store iGUIDE inputs here
project_paths.add_path(
	key = 'raw_sequencing_data', 
	path = pjoin(project_paths.paths['data_raw'], '0-iGUIDE_runs', sequencing_run_id)
	)

# Copy iGUIDE output here
project_paths.add_path(
	key = 'processed_iguide', 
	path = pjoin(project_paths.paths['data_processed'], '0-iGUIDE_runs', sequencing_run_id)
	)

os.makedirs(project_paths.paths['raw_sequencing_data'], exist_ok=True)