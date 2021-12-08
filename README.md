SCATA - Sequence Clustering and Analysis of Tagged Amplicons
============================================================

SCATA is a system for handling and clustering of sequencing data of
tagged amplicon sequences produced through high throughput sequencing
methods. It provides a web frontend where users can upload data
files and select what analyses to perform, as well as a backend that
performs the analysis.


Sequence quality control and demultiplexing
-------------------------------------------

SCATA provides several ways to check quality of the amplicon sequences
before adding them to the analysis. Alla parameters are customisable
thorugh the web interface. If tag sequences are provided, amplicons
are grouped by tag where the tag sequence can be identified.

Sequence clustering
------------------

Sequence clustering in SCATA is performed using highly parallellised
single linkage clustering. Currently two different search enginges can
be used for the clustering process. BLAST and USEARCH. The search
engines are used to identify candidates for inclusion in
clusters. Final scoring of aligments is done by SCATA. One major
feature of the SCATA clustering process is that all sequences of all
samples are clustered simultaneously. Once clusters are identified for
the full experiement with all samples, clusters within each sample are
reconstructed based on global clusters. This workflow ensures that
clusters ("OTUs") can be tracked reliably across samples without any
need of post hoc reconciliation of different per sample clustering
runs.

Furthermore, during the clustering process, reference sequences are
clustered along with all amplicon sequences. The major advantage of
this approach is that clustering settings can be evaluated based on
reference inclusion. For example, over-clustering will need to
multiple references included in the same cluster. More importantly,
reference assigment is done on exactly the same premises as
clustering. Thus, if a reference is included, this is strong evidence
for the conclusion that the cluster and the reference are the same
OTU, given the current clustering settings. However, it is important
to take in to account that references are treated slightly
differently. Most importantly, reference sequences can not form a
bridge between two clusters, that forces a merge of the two clusters
into one. The reference sequence will, rather, be added to both
clusters. Thus, in some rare cases, under-clustering will result in
multiple clusters containing the same reference. 


PRINCIPLE OF OPERATION
======================

The SCATA system consists of a web interface, where users can submit
data and initiate analyses. All user input is stored/queued to an SQL
database, where jobs are set up and made ready to run. The second part
of the system is the backend which is doing the actual heavy lifting,
and dispatches all jobs through a grid middleware (currently SGE). The
backend dispatcher regularly checks database tables for new tagsets,
datasets, reference sets and jobs. When a new entity is found, an SGE
job is launched to check/run the job. The backend dispatcher keeps
track of running jobs, checks the exit status and resubmits the job is
there was a temporary failure. In some cases, the resubmission is
tried with a higher memory allocation request. Tagsets, reference sets
and datasets are all handle as a single job each. A clustering job,
however, is slightly more complex. The main ScataJob will submit many
smaller worker-jobs to cluster the data in parallell, as well as to
merge all subclusterd datasets. Finally, one job for each cluster will
be launched to calculate summary statistics for the job.

INSTALLATION
============

To install SCATA general linux administration experience is assumed,
including setting up PHP and a web server to serve php web
pages. Basic knowledge of how to manage a MariaDB/MySQL database is
also required. SCATA also requires a grid middleware/queue
system. Currently the only supported system is SGE, but it should be
fairly simple to add support for eg SLURM by replacing sge.py with
module written for the specific middleware.

Dependencies
------------

SCATA is written Python3 and uses Biopython to handle sequence
data. There are a number of external dependencies that need to be
installed for a full SCATA web service. The python dependencies are
most conveniently installed within a separate python venv, to ensure
a stable python module environment for SCATA:

 - Python3 (Usually the system python3 works fine)
 - BioPython
 - PyMySQL
 - NCBI BLAST+
 - usearch
 - Web server with "LAMP" stack. 

Prepeare the filesystem
-----------------------

Prepare a filesystem structure like below. You are free to modify it
though, as long as the changes you do are reflected in the
configuration files mention further down in this instruction. 

$ SCATA_ROOT=/scata
$ for i in scata-run/log scata-run/run scata-run/tmp scata-system/bin \
  scata-system scata-files/tmp \
  scata-files/files scata-data/tagsets scata-data/referencesets \
  scata-data/results scata-data/datasets; \
  do mkdir -p $SCATA_ROOT/$i; done

Clone the scata github repository into the scata-system folder.
   
Set up Web interface
--------------------

 1. Prepare a "LAMP" compliant web server.
 2. Set the document root to the www/www directory within the
    distribution.
 3. Configure scata/www/includes/constants.php and set paths according
    to what was created in the previous step.

Set up SQL database
-------------------

 1. Create a database with associated read/write user for SCATA and for
    the SCATA web interface.
 2. Load the scata.schema into the database.
 3. Replace the crypted admin password and username with to your
    preference.
 4. Insert the credentials into the php configuration.

Configure the SCATA dispatcher/daemon
-------------------------------------

The scata system must not be run as root. Please create a separate
scata user, which will own all scata files. The scatad.sh deamon must
be run as this user.

1. Create a Python virtual environment where biopython and PyMySQL are installed.
2. Update scata-bin/constants.py to reflect your system settings.
3. Update paths to activate the virtual environment in sge.py and scatad.sh
4. Set up scatad.sh to start at boot.

