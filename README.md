# nf-core-beav
NF-CORE-BEAV is a bioinformatics pipeline that takes bacteria sequences, either raw fasta sequences or from genbank, and annotates them to highlight bacterial elements. The output are 3 annotated sequence files in different common formats, and a Circos plot with the annotations shown. Tools specific to agrobacterium are also optionally available.

Bacterial Element Annotation reVamped (BEAV) has been revamped again to work in a Nextflow pipeline and integrate into nf-core. Most of its modules are designed to run simultaneously with each other, and on many samples at once(although running it with only 1 sample is also possible).

# Workflow
![Beav Workflow](nf_beav_workflow.png)

# Processes
Each process annotates different sequences of a bacterial genome. A general description is provided here. For information about any specific module, check the ![CITATIONS.md](CITATIONS.md) file.

BAKTA: Standard annotator that runs first on all files.

TIGER2: Integrative and Conjugative Element identification

PhiSpy: Prophage identification

IntegronFinder: Integron identification

MacSyFinder: Secretion Systems identification

DefenseFinder: Defense Systems identification

antiSMASH: Biosynthetic Gene Cluster identification

GapMind: Amino Acid Biosynthesis & Carbon Metabolism Pathway identification

Operon-mapper: Operon prediction

HMMSearch: Dif site identification, T-DNA Borders & Overdrive identification(Agrobacterium Pipeline)

Fuzznuc: pip/tts/nod/hrp box identification, tra/vir box identification(Agrobacterium Pipeline)

Blastn: oriT Region identification

Mobsuite: Plasmid Characterization (nf-core)

fastANI: Biovar & Genomospecies classification(Agrobacterium Pipeline)

k-mers: Ti/Ri Plasmid identification, classification & contigs(Agrobacterium Pipeline)

# Setup
NF-CORE-BEAV can be cloned and used directly from GitHub, or accessed via ![nf-co.re/beav](https://nf-co.re/beav). 
```
git clone https://github.com/weisberglab/nf-core-beav
```
BEAV also uses a database located at ![weisberglab/beavDB](https://github.com/weisberglab/beavDB). This database is **automatically** downloaded on first run, but can be cloned manually:
```
git clone https://github.com/weisberglab/beavDB.git
```
Included in the beavDB repository are custom BAKTA databases for the agrobacterium pipeline and some other necessary patches/files.

The location that other necessary databases will be stored is defined by the parameter --beav_dir, and defaults to ./beavDB. Be sure to replace this with an absolute path(to a location with enough storage, if you need it to download the BAKTA/antiSMASH databases), as the pipeline will not be able to find it if it is run from a different directory each time and will reinstall all databases.

Should certain files be missing from the database directory, the pipeline will automatically reinstall it to beav_dir, overwriting everything in the directory, so please do not store anything in the path specified by --beav_dir besides what is in the weisberglab/beavDB repo and/or the automatically downloaded versions of BAKTA/antiSMASH/TIGER2/TXSS DBs or models.

Open the nextflow.config file in your preferred text editor and create a profile(recommended) or change the defaults for at least the following input parameters:

    input: directory of files, single file, comma-delimited list of files (default: ./inputFiles)
    
    outdir: directory where results are stored (default: ./results)
    
    beav_dir: directory where databases are stored (default: ./beav_dir)
    
    bakta_db: directory where your bakta database is stored, if already downloaded (default: auto downloads to ./beav_dir)
    
    antismash_db: directory where your antiSMASH database is stored, if already downloaded (default: auto downloads to ./beav_dir)
    
    tiger_blast_db: location of a fasta file containing BLAST database information for TIGER to confirm ICE insertions (required, NO DEFAULT)
    
    operon_email: an email for operon mapper results to be sent to (default: none)
    
    agro/agrobacterium: a flag that runs the agrobacterium specific pipeline (default: false)

An example profile has been included in the config file for reference.
```
profiles {
    // Example profile(Use with -profile yourProfileName, edit to your liking or add a new profile with your name below it)
    yourProfileName{
        params{
            bakta_db = './path/to/your/bakta/db'
            outdir   = './my/results/path'
        }

        process {
            // set default options to use for all processes
            executor = 'slurm'

            // specify the name of a specific process to give more resources to
            withName: 'BAKTA_BAKTA' {
                memory  = 12.GB
                cpus    = 2
            }
        }

        conda{
            // configure conda options
            enabled = true
        }

        // if you want to use containers instead, edit them here
        docker.enabled = false
            // these two options do the same thing:
        singularity.enabled = false
        /*singularity{
            enabled = false
        }*/

        env {
            // set default environment variables
            SINGULARITY_CACHEDIR = "$HOME/.singularity/cache"
        }

    }
}
```
Please refer to the ![Nextflow Documentation](https://www.nextflow.io/docs/latest/config.html) for more information on setting up profiles in the config file.

Be sure to include the paths to databases you may already have downloaded(bakta, antiSMASH). Keep in mind that TIGER will not run unless a user supplied BLAST database is provided.

# Usage
After setup has been completed, run the pipeline using the following command:
```
bash
nextflow run nf-core-beav \
  -profile <YOUR_PROFILE> \
  --input <sample.fna,sample2.fna / directory_of_samples> \
  --outdir <OUTDIR>
```
Note that every parameter can be set in the profile, but it is recommended that input and outdir are changed every run for easier file organization.
The following parameters are required to run their respective processes/modules:
  --agrobacterium (runs the agrobacterium pipeline)
  --agro (alias for above)
  --tiger_blast_db <path/to/db.fna> (runs the ICE finder)
  --operon_email <email> (required to run operon mapper because it submits a job to a website, excess usage could clutter the system)

# Output
BEAV produces 6 outputs:
  4 Annotated Genome files in various formats:
    .gbk, .gff, .faa, .ffn
  2 Annotated Figures in two formats using Circos:
    .pdf, .png

Intermediate files are stored in the ./work directory created on run. 
