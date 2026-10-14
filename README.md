<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-beav_logo_dark.png">
    <img alt="nf-core/beav" src="docs/images/nf-core-beav_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/beav)
[![GitHub Actions CI Status](https://github.com/nf-core/beav/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/beav/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/beav/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/beav/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/beav/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/beav)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23beav-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/beav)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/beav**(Bacterial Element Annotation reVamped) is a bioinformatics pipeline that takes bacterial sequences, either raw fasta sequences or from genbank, and annotates them to highlight bacterial elements. The output are 3 annotated sequence files in different common formats, and a Circos plot with the annotations shown. Tools specific to agrobacterium are also optionally available.

![Beav Workflow](docs/images/nf-core-beav_workflow.png)

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

### FROM WINDOWS 10/11
If you are using Windows, please set up a Linux distribution with the following steps. Otherwise, skip to the LINUX/UNIX TERMINAL section.

Step 1: Make sure you have at least 5GB of storage. Open Windows powershell and run:
```
wsl --install
wsl.exe --install Ubuntu
```
Follow the instructions to create a user account and password for the Linux distribution.

Step 2: Start the distribution. Using the ~ indicator will place you in the home directory:
```
wsl ~
```
Step 3: Install miniconda3 via the following commands:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Follow its prompts to complete setup.

Step 4: Proceed to the next section to set up Beav.

### FROM LINUX/UNIX TERMINAL
Step 1: Clone the pipeline from [weisberglab/nf-core-beav](https://github.com/weisberglab/nf-core-beav). In future versions, nf-core-beav will be integrated with nf-core and may be called directly from the command line.

Step 2: Set nextflow and environment variables. Open the 'nextflow.config' file in your preferred editor and create a profile(recommended) or change the defaults for the input parameters. 
  Important parameters to note: 
  
    input: directory of files, single file, comma-delimited list of files (default: ./inputFiles)
    outdir: directory where results are stored (default: ./results)
    beav_dir: directory where databases are stored (default: ./beav_dir)
    bakta_db: directory where your bakta database is stored, if already downloaded (default: auto downloads to ./beav_dir)
    antismash_db: directory where your antiSMASH database is stored, if already downloaded (default: auto downloads to ./beav_dir)
    tiger_blast_db: location of a fasta file containing BLAST database information for TIGER to confirm ICE insertions (required, NO DEFAULT)
    operon_email: an email for operon mapper results to be sent to (default: none)
    agro/agrobacterium: a flag that runs the agrobacterium specific pipeline (default: false)

More information on setting up Nextflow can be found in the link contained in the note at the beginning of this Usage section.

Step 3: Install dependencies. It is recommended to use a conda environment. The pipeline SHOULD be able to use either conda or singularity/apptainer to run its modules, so choose one or the other. Execute this command inside the directory containing the pipeline(where the beav_env.yml file is stored) if you choose to run the pipeline with Conda. 
 ```bash
 cd /path/to/nf-core-beav
 conda env create -n beav_env -f beav_env.yml
 ```

Step 4: Run the pipeline. It will automatically set up and download any databases you did not include in the parameters. Note that --input may be a PATH to a single file, a directory containing files(depth 1, automatically ignores unrelated files), or a comma-delimited list of files.
  ```bash
  nextflow run nf-core/beav \
    -profile <YOUR_PROFILE> \
    --input <samplesheet.csv || sample.fna,sample2.fna || directory_of_samples> \
    --outdir <OUTDIR>
  ```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/beav/usage) and the [parameter documentation](https://nf-co.re/beav/parameters).

## Pipeline output


To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/beav/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/beav/output).

## Features

Main Pipeline

BAKTA: Automatic download of Bakta database (light/full) and Bakta annotation, and automatically skip Bakta on input files detected to already be annotated by Bakta

ANTISMASH: Automatic download of AntiSMASH database and AntiSMASH annotation

MOBSUITE_RECON: Reconstruct plasmids in bacterial assemblies

MACSYFINDER: Annotation of macromolecular systems, genetic pathways

DEFENSEFINDER: Annotation of defense systems

GAPMIND: Annotation of amino acid and carbon metabolism

NHMMER: Annotation of structural RNAs, transposable elements, other

FUZZNUC: Annotation of pip/tts/nod box

BLASTN: Search and annotation of oriT, dif sites

SOURMASH: Comparison to reference metagenomes

TIGER: Annotation of ICEs

INTEGRONFINDER: Annotation of integrons

PHISPY: Annotation of prophages

CIRCOS: Create a circular plot with all annotated elements


Agro-specific Pipeline (Optional)

FASTANI: Whole genome ANI computation

FUZZNUC: Annotation of tra/vir box

NHMMER: Annotation of T-DNA borders and Overdrive

IDENTIFY_TI_RI_PLASMIDS: Identify and classify Ti and Ri plasmids and contigs

## Credits

nf-core/beav was originally written by @amnalert, @jewell-bug, @acarafat, and @alexweisberg.

We thank the following people for their extensive assistance in the development of this pipeline:

  Alexandra Weisberg
  Arafat Rahman
  Jewell Jung
  The Oregon State University Center for Qualitative Life Sciences Team
  The Nextflow and nf-core community

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#beav` channel](https://nfcore.slack.com/channels/beav) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/beav for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
