
import os
from os.path import join as pjoin
from src.ProjectPaths import ProjectPaths

# Base path
project_paths = ProjectPaths.ProjectPaths('/data/GenomicTrackRepository')


os.makedirs(pjoin(project_paths.paths['notebooks'], 'hg38'), exist_ok = True)

os.makedirs(pjoin(project_paths.paths['data_processed'], 'hg38'), exist_ok = True)
