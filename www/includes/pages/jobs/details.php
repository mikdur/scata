<h1>Job Details</h1>
<table>
<?php
if (isset($_GET['id']) && $_GET['id'] != "")
{
	$job = new Job($_GET['id']);
	if (!$job->user_has_rights())
	{
		redirect("./?p=jobs");
	}

	// Beskrivning
	echo "<tr><td><b>Name:</b></td><td>".$job->get_name();
	echo "</td></tr>";

	// Beskrivning
	echo "<tr><td><b>Description:</b></td><td>".$job->get_description();
	echo "</td></tr>";
	
	// Skriv ut parameterset
	$paramset = new ParameterSet($job->get_parameterset());
	echo "<tr><td><b>Parameter set: </b></td><td>".$paramset->get_name();
	echo "</td></tr>";
	

	// Skriv ut dataset
	$datasets = $job->get_datasets();
	echo "<tr><td><b>Datasets:</b></td><td>";
	echo "<ul>";
	foreach($datasets as $dataset)
	{
		$dataset = new Dataset($dataset);
		echo "<li>".$dataset->get_name()."</li>";
	}
	echo "</ul></td></tr>";
	
	// Skriv ut referenssekvenser
	$refseqs = $job->get_refseqs();
	echo "<tr><td><b>Reference sequences:</b></td><td>";
	echo "<ul>";
	foreach($refseqs as $refseq)
	{
		$refseq = new Refseq($refseq);
		echo "<li>".$refseq->get_name()."</li>";
	}
	echo "</ul></td></tr>";

}
else
{
	redirect("./?p=jobs");
}
?>
</table>
<a href="./?p=jobs">Back</a>