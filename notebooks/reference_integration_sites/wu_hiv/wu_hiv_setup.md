Alex McFarland

10/22/2024

# Background

```sh
PROJECT_PATH="/data/GenomicTrackRepository"
cd $PROJECT_PATH
conda activate repeat_project
ipython
```

# Import libraries

```python
import os
import pandas as pd
from os.path import join as pjoin
import glob
import gzip
from Bio import SeqIO
from src.ProjectPaths import ProjectPaths

pd.set_option('display.max_columns', None)

project_paths = ProjectPaths.ProjectPaths('/data/GenomicTrackRepository')

analysis_name = 'wu_hiv'
```

# Set working paths

```python
project_paths.add_path(
	key = 'data_raw_sequencing', 
	path = pjoin(project_paths.paths['data_raw'], 'reference_integration_sites', analysis_name)
	)
os.makedirs(project_paths.paths['data_raw_sequencing'], exist_ok = True)

project_paths.add_path(
	key = 'analysis_interim', 
	path = pjoin(project_paths.paths['data_interim'], 'reference_integration_sites', analysis_name)
	)
os.makedirs(project_paths.paths['analysis_interim'], exist_ok = True)


project_paths.add_path(
	key = 'anaylsis_processed', 
	path = pjoin(project_paths.paths['data_processed'], 'reference_integration_sites', analysis_name)
	)
os.makedirs(project_paths.paths['anaylsis_processed'], exist_ok = True)

aavenger_dir = '/data/AAVengeR'
```

# Add metadata to raw folder

```python

df_metadata = pd.DataFrame({
    "run_date": [
        "2018-03-05", "2019-02-01", "2019-09-12", "2019-09-15", "2019-09-20", 
        "2019-12-10", "2019-12-18", "2019-12-20"
    ],
    "run_ID": [
        "180305_M00281_0327_000000000-BMK22", "190201_M00281_0451_000000000-C8N97", 
        "190912_M03249_0010_000000000-CKNGK", "190915_M00281_0533_000000000-CN99V", 
        "190920_M03249_0013_000000000-CNJ7H", "191210_M05588_0259_000000000-CRV2M", 
        "191218_M03249_0040_000000000-CRT4T6", "191220_M03249_0042_000000000-CRRWT"
    ],
    "instrument": ["Hahn", "Hahn", "Bushman", "Hahn", "Bushman", "Hahn", "Bushman", "Bushman"],
    "fastq_path_microb120": [
        "/home/wuv/int1_data/fastq", "/home/wuv/int5_data/fastq", 
        "/home/wuv/int6_data/A-C_fastq", "/home/wuv/int6_data/D-F_fastq", 
        "/home/wuv/int6_data/G-I_fastq", "/home/wuv/int7_data/A-C", 
        "/home/wuv/int7_data/D-F", "/home/wuv/int7_data/G-I"
    ],
    "sample_info": [
        "ART-treated; run 1", "ART-treated; run 2", 
        "acute infection; samples A-C", "acute infection; samples D-F", 
        "acute infection; samples G-I", "chronic infection; samples A-C", 
        "chronic infection; samples D-F", "chronic infection; samples G-I"
    ],
    "int_site_exp": [1, 5, 6, 6, 6, 7, 7, 7]
})

df_metadata['raw_data_local'] = df_metadata.apply(lambda x: pjoin(project_paths.paths['data_raw_sequencing'], x['run_ID']), axis = 1)


df_metadata['interim_data_local'] = df_metadata.apply(lambda x: pjoin(project_paths.paths['analysis_interim'], x['run_ID']), axis = 1)

df_metadata['processed_data_local'] = df_metadata.apply(lambda x: pjoin(project_paths.paths['anaylsis_processed'], x['run_ID']), axis = 1)

print(df_metadata)
```


# Copy data from microb120 to aws

Copying the sample sheet from Vincent Wu's project directory. At the moment no int6_data csv metadata is in his folder.

Manually populated folders for int6_data from emails on 10/22/24

```python
for _, run_ in df_metadata.iterrows():

	os.makedirs(run_.raw_data_local, exist_ok = True)

	os.system(f"scp agmcfarland@microb120.med.upenn.edu:{run_.fastq_path_microb120}/*fastq.gz {run_.raw_data_local}")

	os.system(f"scp agmcfarland@microb120.med.upenn.edu:{run_.fastq_path_microb120}/*.csv {run_.raw_data_local}")
```


# Gather run metadata

```python
sample_sheet_supplemental = []

sample_sheet_raw_path = []

for _, run_ in df_metadata.iterrows():

	for f in glob.glob(pjoin(run_.raw_data_local, '*.csv')):
		if f.endswith('.supp.csv'):
			sample_sheet_supplemental.append(f)
		else:
			sample_sheet_raw_path.append(f)

df_metadata['sample_sheet_raw_path'] = 	sample_sheet_raw_path

df_metadata['sample_sheet_supplemental_path'] = sample_sheet_supplemental
```

# Check how well hmm works

HIV-1 PCR2 (U3), 

CAAGCAGAAGACGGCATACGAGATBARCODEAGTCAGTCAGCCCAGGGAAGTAGC CTTGTGTGTGGT

CTTGTGTG -- primer tip

HIV-1 PCR2 (U5), 

CAAGCAGAAGACGGCATACGAGATBARCODEAGTCAGTCAGCCCAAGTAGTGTGTGCCC GTCTGTTG

GTCTGTTG -- primer tip

Using default AAVengeR HMMs and settings for wild HIV

```python

for _, run_ in df_metadata.iterrows():break

	with gzip.open(pjoin(run_.raw_data_local, 'Undetermined_S0_R2_001.fastq.gz'), 'rt') as input_handle:
		store_record = []
		for record in SeqIO.parse(input_handle, 'fastq'):
			store_record.append(record)

	for f in glob.glob(pjoin(run_.raw_data_local, '*.csv')):
		if f.endswith('.supp.csv'):
			sample_sheet_supplemental.append(f)
		else:
			sample_sheet_raw_path.append(f)

df_metadata['sample_sheet_raw_path'] = 	sample_sheet_raw_path

df_metadata['sample_sheet_supplemental_path'] = sample_sheet_supplemental


u3_records = []
for record in store_record:
	if str(record.seq).startswith('CTTGTGTGTGGT'):
		u3_records.append(record)

u3_first_x_bp = []
for record in u3_records:
	u3_first_x_bp.append(str(record.seq)[:100])


for record in u3_records:
	if str(record.seq).find('CCCTTCCA')


with open(pjoin(project_paths.paths['analysis_interim'], 'U3_subset.fasta'), 'wt') as outfile:
	for record in u3_records:
		outfile.write(f'>{record.description}\n{str(record.seq)}\n')

os.system(f"nhmmer --cpu 8 --max --tblout {pjoin(project_paths.paths['analysis_interim'], 'U3_subset.tsv')} {pjoin(aavenger_dir, 'data/hmms', 'HIV1_1-100_U3_RC.hmm')} {pjoin(project_paths.paths['analysis_interim'], 'U3_subset.fasta')} > {pjoin(project_paths.paths['analysis_interim'], 'U3_subset.stdout')}")


u5_records = []
for record in store_record:
	if str(record.seq).startswith('GTCTGTTG'):
		u5_records.append(record)

u5_first_x_bp = []
for record in u5_records:
	u5_first_x_bp.append(str(record.seq)[:100])

with open(pjoin(project_paths.paths['analysis_interim'], 'u5_subset.fasta'), 'wt') as outfile:
	for record in u5_records:
		outfile.write(f'>{record.description}\n{str(record.seq)}\n')

os.system(f"nhmmer --cpu 8 --max --tblout {pjoin(project_paths.paths['analysis_interim'], 'u5_subset.tsv')} {pjoin(aavenger_dir, 'data/hmms', 'HIV1_1-100_U5.hmm')} {pjoin(project_paths.paths['analysis_interim'], 'u5_subset.fasta')} > {pjoin(project_paths.paths['analysis_interim'], 'u5_subset.stdout')}")

```

# Make complete metadata sheet for aavenger

```python
for _, run_ in df_metadata.iterrows():

	sample_info = pd.read_csv(run_.sample_sheet_supplemental_path)

	sample_info = sample_info.rename(columns = {'specimen': 'subject'})
	
	sample_sheet_raw = pd.read_csv(run_.sample_sheet_raw_path)

	sample_sheet_raw['sample'] = sample_sheet_raw['sampleName']

	# remove uninfected and no template controls
	sample_sheet_raw = sample_sheet_raw[~sample_sheet_raw['sample'].str.contains('infected|emplate|UNC|NTC')]

	sample_sheet_raw['subject_1'] = sample_sheet_raw['sample'].apply(lambda x: x.split('-')[0])	

	# Group by the two columns and assign IDs 
	sample_sheet_raw['replicate'] = sample_sheet_raw.groupby(['subject_1', 'uniqueRegion']).cumcount() + 1 #% 4 + 1

	sample_sheet_raw['trial'] = run_.run_ID

	sample_sheet_raw['vectorFastaFile'] = 'HXB2.fasta'
	
	sample_sheet_raw['flags'] = sample_sheet_raw['uniqueRegion'].apply(lambda x: 'IN_u5' if x.find('U5') > -1 else 'IN_u3')
	
	sample_sheet_raw['leaderSeqHMM'] = sample_sheet_raw['uniqueRegion'].apply(lambda x: 'HIV1_1-100_U5.hmm' if x.find('U5') > -1 else 'HIV1_1-100_U3_RC.hmm')
	
	sample_sheet_raw['adriftReadLinkerSeq'] = sample_sheet_raw['barcode1']

	sample_sheet_raw['index1Seq'] = sample_sheet_raw['barcode2']

	sample_sheet_raw['refGenome'] = 'hs1'

	# Add metadata

	shape_prior_to_merge = sample_sheet_raw.shape

	sample_sheet_completed = sample_sheet_raw.merge(
		sample_info,
		left_on = 'subject_1',
		right_on = 'subject',
		how = 'outer'
		)

	# expect no samples weren't lost or gained when merging metadata
	assert shape_prior_to_merge[0] == sample_sheet_completed.shape[0]

	# subset just columns used by AAVengeR

	sample_sheet_completed = sample_sheet_completed[['trial', 'subject', 'sample', 'replicate', 'adriftReadLinkerSeq', 'index1Seq', 'refGenome', 'vectorFastaFile', 'leaderSeqHMM', 'flags']]

	# expect each sample to have four replicates each for U3 and U5
	assert(len(sample_sheet_completed)/4//2 == len(sample_sheet_completed['subject'].unique()))

	# expect all sample names to be unique

	assert(len(sample_sheet_completed['sample'].unique()) == len(sample_sheet_completed['sample'].tolist()))

	sample_sheet_completed.to_csv(pjoin(run_.raw_data_local, 'CompletedSampleSheet.tsv'), sep = '\t', index = None)
```


# Make config file

Manually added config file to external data. Added time stamp for reference.

```python
production_config = pjoin(project_paths.paths['data_external'], 'production_config_20241022.yml')

for _, run_ in df_metadata.iterrows():

	config_file_completed = pjoin(run_.raw_data_local, 'CompletedConfigFile.yml')

	with open(config_file_completed, 'w') as outfile:

		with open(production_config, 'r') as infile:

			for l in infile:

				if l.startswith('outputDir'):
					outfile.write(f"outputDir: {run_.processed_data_local}\n")

				elif l.startswith('sequencingRunID'):
					outfile.write(f"sequencingRunID: {run_.run_ID}\n")

				elif l.startswith('demultiplex_anchorReadsFile'):
					outfile.write(f"demultiplex_anchorReadsFile: {run_.raw_data_local}/Undetermined_S0_R2_001.fastq.gz\n")

				elif l.startswith('demultiplex_adriftReadsFile'):
					outfile.write(f"demultiplex_adriftReadsFile: {run_.raw_data_local}/Undetermined_S0_R1_001.fastq.gz\n")

				elif l.startswith('demultiplex_index1ReadsFile'):
					outfile.write(f"demultiplex_index1ReadsFile: {run_.raw_data_local}/Undetermined_S0_I1_001.fastq.gz\n")

				elif l.startswith('demultiplex_sampleDataFile'):
					outfile.write(f"demultiplex_sampleDataFile: {run_.raw_data_local}/CompletedSampleSheet.tsv\n")

				elif l.startswith('annotateRepeats_inputFile: '):
					outfile.write('annotateRepeats_inputFile: core/sites.rds\n')

				elif l.startswith('callNearestGenes_inputFile: '):
					outfile.write('callNearestGenes_inputFile: annotateRepeats/sites.rds\n')

				elif l.find('_CPUs: ') > -1:
					l = l.replace('_CPUs: ', '_CPUs: 30 #')
					outfile.write(f"{l}\n")

				elif l.startswith('softwareDir: '):
					outfile.write('softwareDir: /data/AAVengeR\n')

				else:
					outfile.write(l)
```


# Run AAVengeR with docker

```sh
screen -S aavenger_run

sudo docker run -it --mount type=bind,source=/data,dst=/data/ aavenger_docker_v2_2 bash 

cd /data/GenomicTrackRepository/data/raw/reference_integration_sites/wu_hiv

AAVENGER_DIR="/data/AAVengeR"

for config_ in /data/GenomicTrackRepository/data/raw/reference_integration_sites/wu_hiv/*/CompletedConfigFile.yml
do
	# head -n 100  $config_
	Rscript $AAVENGER_DIR/aavenger.R $config_
done


```

```sh
# monitoring
cd /data/GenomicTrackRepository/data/processed/reference_integration_sites/wu_hiv/180305_M00281_0327_000000000-BMK22/core
tail -f demultiplex/log
tail -f replicateJobTable
```

# Remove raw sequencing files

```python
[os.remove(f) for f in  glob.glob(pjoin(project_paths.paths['data_raw_sequencing'], '*fastq.gz'))]
````

df_run_metadata = pd.read_csv(pjoin(project_paths.paths['data_processed'], '0_process_aavenger_runs', 'aavenger_run_metadata.csv'))

df_run_metadata = df_run_metadata[df_run_metadata['run_id'] == illumina_run_id]

# emi_study_id = '_'.join(list(df_run_metadata['Emi_StudyID'].unique()))
emi_study_id = 'limit_test'

print(df_run_metadata['GTSP'].unique())

print(emi_study_id)
```

# Pull data

```python
os.system(f"scp agmcfarland@microb120.med.upenn.edu:/media/sequencing/Illumina/{illumina_run_id}/Data/Intensities/BaseCalls/*.fastq.gz {project_paths.paths['data_raw_sequencing']}")

os.system(f"scp agmcfarland@microb120.med.upenn.edu:/media/sequencing/Illumina/{illumina_run_id}/SampleSheet.csv {project_paths.paths['data_raw_sequencing']}")
```

# Make sample sheet for AAVengeR

```python
sample_sheet_raw = pd.read_csv(pjoin(project_paths.paths['data_raw_sequencing'], 'SampleSheet.csv'), skiprows = 23, usecols = [0, 1, 2, 3])

sample_sheet_raw['sample'] = sample_sheet_raw['SampleName'].apply(lambda x: x.split('-')[0])
sample_sheet_raw = sample_sheet_raw.sort_values(['sample', 'uniqueRegion'], ascending = [True, False])

sample_sheet_raw['replicate'] = [1, 2, 3, 4, 5, 6, 7, 8] * (sample_sheet_raw.shape[0]//8)
sample_sheet_raw['trial'] = emi_study_id
sample_sheet_raw['vectorFastaFile'] = 'HXB2.fasta'
sample_sheet_raw['flags'] = sample_sheet_raw['SampleName'].apply(lambda x: 'IN_u5' if x.find('_u5') > -1 else 'IN_u3')
sample_sheet_raw['leaderSeqHMM'] = sample_sheet_raw['SampleName'].apply(lambda x: 'HIV1_1-100_U5.hmm' if x.find('_u5') > -1 else 'HIV1_1-100_U3_RC.hmm')
sample_sheet_raw['adriftReadLinkerSeq'] = sample_sheet_raw['linker']
sample_sheet_raw['index1Seq'] = sample_sheet_raw['barcode']
sample_sheet_raw['refGenome'] = 'hg38'

sample_sheet_raw = sample_sheet_raw.merge(
	df_run_metadata[['full_GTSP', 'Bushman_MouseID']],
	left_on = 'SampleName',
	right_on = 'full_GTSP',
	how = 'inner')

sample_sheet_raw['subject'] = sample_sheet_raw['Bushman_MouseID']

print(len(sample_sheet_raw['Bushman_MouseID'].unique()))

assert(len(sample_sheet_raw)/4//2 == len(sample_sheet_raw['sample'].unique()))

sample_sheet_completed = sample_sheet_raw[['trial', 'subject', 'sample', 'replicate', 'adriftReadLinkerSeq', 'index1Seq', 'refGenome', 'vectorFastaFile', 'leaderSeqHMM', 'flags']]

sample_sheet_completed.to_csv(pjoin(project_paths.paths['data_raw_sequencing'], 'CompletedSampleSheet.tsv'), sep = '\t', index = None)
```

# Make config file

```python
production_config = pjoin(project_paths.paths['data_external'], 'production_config.yml')

config_file_completed = pjoin(project_paths.paths['data_raw_sequencing'], 'CompletedConfigFile.yml')

with open(config_file_completed, 'w') as outfile:

	with open(production_config, 'r') as infile:

		for l in infile:

			if l.startswith('outputDir'):
				outfile.write(f"outputDir: {project_paths.paths['aavenger_results']}\n")

			elif l.startswith('demultiplex_seqRunID'):
				outfile.write(f"demultiplex_seqRunID: {illumina_run_id}\n")

			elif l.startswith('demultiplex_anchorReadsFile'):
				outfile.write(f"demultiplex_anchorReadsFile: {project_paths.paths['data_raw_sequencing']}/Undetermined_S0_R2_001.fastq.gz\n")

			elif l.startswith('demultiplex_adriftReadsFile'):
				outfile.write(f"demultiplex_adriftReadsFile: {project_paths.paths['data_raw_sequencing']}/Undetermined_S0_R1_001.fastq.gz\n")

			elif l.startswith('demultiplex_index1ReadsFile'):
				outfile.write(f"demultiplex_index1ReadsFile: {project_paths.paths['data_raw_sequencing']}/Undetermined_S0_I1_001.fastq.gz\n")

			elif l.startswith('demultiplex_sampleDataFile'):
				outfile.write(f"demultiplex_sampleDataFile: {project_paths.paths['data_raw_sequencing']}/CompletedSampleSheet.tsv\n")

			elif l.startswith('annotateRepeats_inputFile: '):
				outfile.write('annotateRepeats_inputFile: core/sites.rds\n')

			elif l.startswith('callNearestGenes_inputFile: '):
				outfile.write('callNearestGenes_inputFile: annotateRepeats/sites.rds\n')

			elif l.find('_CPUs: ') > -1:
				l = l.replace('_CPUs: ', '_CPUs: 30 #')
				outfile.write(f"{l}\n")

			elif l.startswith('softwareDir: '):
				outfile.write('softwareDir: /data/AAVengeR\n')

			else:
				outfile.write(l)
```


# Run AAVengeR with docker

```sh
screen -S aavenger_run

sudo docker run -it --mount type=bind,source=/data,dst=/data/ aavenger_docker_v2_2 bash 

AAVENGER_DIR="/data/AAVengeR"

CONFIG_PATH="/data/brad_jones_cd8_pdx_project/cd8_pdx_project/data/raw/aavenger/240912_MN01490_0259_A000H7GFVK/CompletedConfigFile.yml"
 
Rscript $AAVENGER_DIR/aavenger.R $CONFIG_PATH
```

```sh
#monitoring
cd /data/brad_jones_cd8_pdx_project/cd8_pdx_project/data/interim/aavenger/240912_MN01490_0259_A000H7GFVK/core
tail -f demultiplex/log
tail -f replicateJobTable
```

# Remove raw sequencing files

```python
[os.remove(f) for f in  glob.glob(pjoin(project_paths.paths['data_raw_sequencing'], '*fastq.gz'))]
````

# Check quality

```sh
conda activate vittoria_env

cd /data/brad_jones_cd8_pdx_project/cd8_pdx_project/data/raw/aavenger/240912_MN01490_0259_A000H7GFVK

fastqc Undetermined_S0_R1_001.fastq.gz Undetermined_S0_R2_001.fastq.gz Undetermined_S0_I1_001.fastq.gz

rm *_fastqc.zip
```









