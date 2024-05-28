
import os
from os.path import join as pjoin
from src.ProjectPaths import ProjectPaths

# Base path
project_paths = ProjectPaths.ProjectPaths('/data/GenomicTrackRepository')

for genome_assembly_ in ['hg38', 'mm9']:

	os.makedirs(pjoin(project_paths.paths['notebooks'], genome_assembly_), exist_ok = True)

	os.makedirs(pjoin(project_paths.paths['data_processed'], genome_assembly_), exist_ok = True)
