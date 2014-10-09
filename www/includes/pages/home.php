<h1>SCATA - Sequence Clustering and Analysis of Tagged Amplicons</h1>

<p>
SCATA provids an analysis framework for the analysis of sequenced tagged 
amplicons, typically derived from high throughput sequencing of
microbial communities. It is optimised for target sequences which cannot
readily be aligned across wide phylogenies, e.g. the ITS region. For multiple alignable
target sequences, such as 16S rRNA, we recommend the use of purposebuilt
systems, eg. tools provided by the <a href="http://rdb.cme.msu.edu">Ribosomal Database Project</a>.
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
<b>Wiki with documentation</b><br>
Scata documentation is growing within the <a href="https://scata.mykopat.slu.se/trac">Scata Wiki</a>. 
Not all documentation is up to date, but we are working on it.
If you find Scata a useful service, and use it in a way that can be useful for other people, please feel
free to contribute your usage case to the wiki. (2011-06-20)
</li>
<li>
<b>SCATA updated!</b><br/>We have updated several aspects of SCATA during the last few
months. These changes affects several aspects of data handling under the hood
as well as a number of things affecting you as users. Please se below for some
details of what has changed. The old version of SCATA is no longer available (2011-03-23)
</li>
<li>
<b>Improved dataset upload</b><br/>
We have moved the quality screening and filtering of datasets to the import screen.
This implies that you only have to wait for this process once per dataset, instead of
every time you run the clustering. Due to this, all dataset which have previously
been uploaded have to be uploaded again. Another new feature for the filtering
is the ability to extract long high quality regions from the reads, instead of
basing the filtering on the complete read. This can be good for datasets which are not trimmed
by the Roche software. SCATA is now also able to detect primer sequences
in either end of amplicons in order reverse complement sequences as necessary. We have 
also added support for tags in both ends. (2011-03-23)
</li>
<li>
<b>Change of homopolymer handling</b><br/>
Scata now gives the option to collapse
homopolymers over a given length before clustering. This only affects the sequences
in the search and clustering process. Sequences in the report files are not 
homopolymer collapsed and thus fully comparable to references in other databases. (2011-03-23)
</li>
<li>
<b>Change of default parameters</b><br/>
Many parameters in parameter sets have got
new default values. Most notably, with the new homopolymer handling, we recommend the
use of gap open and extension penalties to avoid over clustering. (2011-03-23)
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
