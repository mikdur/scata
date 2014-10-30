<h1>SCATA - Sequence Clustering and Analysis of Tagged Amplicons</h1>

<p>
SCATA provids an analysis framework for the analysis of sequenced tagged 
amplicons, typically derived from high throughput sequencing of
microbial communities. It is optimised for target sequences which cannot
readily be aligned across wide phylogenies, e.g. the ITS region. For multiple alignable
target sequences, such as 16S rRNA, we recommend the use of pipelines
optimised for such data.
</p>

<p>
Please note that the Scata service is offered freely to the non-commercial 
scientific community, and as such is run on otherwise unused computer
time. This implies that at times, analyses will take longer (up to several days)
to finish depending on other requirments of other projects for computational resources. 
</p>
<h2>News</h2>
<p>
<ul>
  <li>
    <b>Major update of scata capabilites!</b>
    Today we launched a major revision of the scata system with
    support for recent sequencing advances. This includes support for
    reading files from both IonTorrent and MiSeq systems. More
    specifically the following things have been added/updated:
    <ul>
      <li><b>New file format support</b> - Scata now supports uploading
      of sanger fastq files, in addition to the previously supported
	formats (fasta, fasta+qual, roche/454 sff). The fastq files
      can be gzipped, which decreases the upload time.</li>
      <li><b>Merging of paired MiSeq reads</b> - Scata can now merge
      paired MiSeq reads by aligning the reads and checking for
	overlap. </li>
      <li><b>Amplicon quality filtering</b> - New filter, which checks
      the quality of the amplicon once the primers have been
      identified. Useful for sequencing methods where the quality
      drops towards the end, and where low quality bases can occur
	within the sequence.</li>
      <li><b>Separate upload and dataset verification</b> - Datasets
      should now be uploaded under the File tab, and can then be
      checked and imported as dataset. This makes it easier to try
      different quality check parameters without having to send the
	whole file again.</li>
      <li><b>Dataset check speedup</b> - Checking of datasets is now
	5-10 times faster than before!</li>
      <li><b>Full support for ambiguity bases in primer sequences</b>
      - Scata now has full support for the IUPAC ambiguity codes in
      the primer sequence. Thus, you should probably leave the primer
      score at 0.9 or similar even in cases where there are ambiguity
	bases in the primer sequence when importing new datasets.</li>
      <li><b>Fixed 3´ tag support</b> - There was a bug in the 3´
	primer matching code which is now fixed.</li>
      <li><b>Bug fixes</b> - There are several other minor bugfixes
      and cleanups to the code which improves stability, speed and
      ease of use.</li>
      </ul>
    (2014-10-28)
  </li>
<li>
<b>Wiki with documentation</b><br>
Scata documentation is growing within the <a href="https://scata.mykopat.slu.se/trac">Scata Wiki</a>. 
Not all documentation is up to date, but we are working on it.
If you find Scata a useful service, and use it in a way that can be useful for other people, please feel
free to contribute your usage case to the wiki. (2011-06-20)
</li>

</ul>

</p>
<h2>Disclaimer</h2>
<p>We provide SCATA as a free service to the scientific community. Please
make sure to download your results when done, as we cannot take any long-term
responsibility for data storage (datasets are large and use hard disk real 
estate!). We have tested the analysis pipeline throroughly and use it 
regularly for our own projects. However, we cannot guarantee that it is error-free; 
the final responsibility for ensuring that your results are correct
rest with you.</p>

<p>We make no warranty (expressed, implied, or statutory) regarding any
data stored whithin this service or any results obtained through using it,
including without limitation implied warranties of merchantability, 
fitness for use, or fitness for a particular purpose.</p>

 

<h2>Citing SCATA</h2>
<p>
A paper describing SCATA is now submitted. The current citation
if you use this service is:<br/>
<pre>
Mikael Brandström Durling, Karina E Clemmensen, Jan Stenlid and Björn
Lindahl (2011): SCATA - An efficient bioinformatic pipeline for species identification and quantification after high-throughput sequencing of tagged amplicons (submitted).
</pre>
</p>
<p>
This service is sponsored by the <a href="http://www.mykopat.slu.se/">
Department of Forest Mycology and Plant Pathology</a> at the <a
href="http://www.slu.se">Swedish University of Agricultural Sciences.</a>  
Please direct any questions regarding the system to Mikael Brandström Durling
using the email address mikael::durling@slu:;se (replacing the colons 
with dots).
</p>




<?php  #phpinfo() ?>
